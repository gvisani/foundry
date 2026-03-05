# '''
# Tailored dataset wrappers for design tasks
# '''

import json
import os
import textwrap
from os import PathLike
from typing import Any, Dict, List

import numpy as np
import torch
import yaml
from atomworks.ml.datasets import MolecularDataset
from atomworks.ml.transforms.base import Compose, Transform
from omegaconf import DictConfig, OmegaConf
from rfd3.inference.input_parsing import (
    DesignInputSpecification,
    ensure_input_is_abspath,
)
from rfd3.inference.parsing import InputSelection
from torch.utils.data import (
    DataLoader,
    SequentialSampler,
)

from foundry.utils.alignment import weighted_rigid_align
from foundry.utils.datasets import assemble_distributed_loader
from foundry.utils.ddp import RankedLogger

logger = RankedLogger(__name__, rank_zero_only=True)
all_ranks_logger = RankedLogger(__name__, rank_zero_only=False)


def _align_alt_to_main(data, data_alt, align_on):
    """
    Align alt coordinates to main coordinates using rigid body alignment.
    Modifies data_alt's atom_array coordinates in-place.

    Must be called BEFORE transforms, so that residue numbering is preserved
    and InputSelection can resolve residue IDs from the original PDB.

    Args:
        data: main pre-transform pipeline output dict (has "atom_array")
        data_alt: alt pre-transform pipeline output dict (has "atom_array")
        align_on:
            True  → align on all shared fixed CA atoms matched by (chain_id, res_id)
            dict  → e.g. {"A1-100,B1-50": "C40-140,D1-50"} maps main selections to alt selections;
                    only CA atoms within each selection are used for alignment
    """
    aa_main = data["atom_array"]
    aa_alt = data_alt["atom_array"]


    if isinstance(align_on, dict):
        # --- Dict mapping mode ---
        # Each key selects atoms from main, each value selects atoms from alt.
        # Only CA atoms within each selection are used. Counts must match.
        main_indices = []
        alt_indices = []

        for main_sel_str, alt_sel_str in align_on.items():
            main_sele = InputSelection.from_any(main_sel_str, atom_array=aa_main)
            alt_sele = InputSelection.from_any(alt_sel_str, atom_array=aa_alt)

            main_idx = np.where(main_sele.mask & (aa_main.atom_name == "CA"))[0]
            alt_idx = np.where(alt_sele.mask & (aa_alt.atom_name == "CA"))[0]

            assert len(main_idx) == len(alt_idx), (
                f"Alignment selection mismatch: main '{main_sel_str}' resolved to "
                f"{len(main_idx)} CA atoms, alt '{alt_sel_str}' resolved to {len(alt_idx)} CA atoms. "
                f"Both sides must select the same number of CA atoms."
            )

            main_indices.extend(main_idx.tolist())
            alt_indices.extend(alt_idx.tolist())

        assert len(main_indices) > 0, "No CA atoms selected for alignment."

    else:
        raise ValueError(f"Invalid align_on value: {align_on!r}. Must be a dict mapping main selections to alt selections.")

    # --- Compute rigid alignment transform from matched atoms ---
    coords_main = torch.tensor(aa_main.coord[main_indices], dtype=torch.float32)  # [N, 3]
    coords_alt = torch.tensor(aa_alt.coord[alt_indices], dtype=torch.float32)     # [N, 3]

    # Compute centers of mass for matched atoms
    center_main = coords_main.mean(dim=0)   # [3]
    center_alt = coords_alt.mean(dim=0)     # [3]

    # Compute rotation via SVD (same method as weighted_rigid_align)
    centered_alt = coords_alt - center_alt    # [N, 3]
    centered_main = coords_main - center_main  # [N, 3]

    C = centered_alt.T @ centered_main  # [3, 3]
    U, S, V = torch.linalg.svd(C)

    # Ensure proper rotation (det = +1)
    F_det = torch.eye(3)
    F_det[-1, -1] = torch.sign(torch.linalg.det(U @ V))
    R = U @ F_det @ V  # [3, 3]

    # Apply transform to ALL alt atom coordinates (modifies atom array in-place)
    all_alt_coords = torch.tensor(aa_alt.coord, dtype=torch.float32)  # [L_alt, 3]
    aligned_all = (all_alt_coords - center_alt) @ R + center_main
    aa_alt.coord = aligned_all.numpy()

    logger.info(
        f"Aligned alt structure to main using {len(main_indices)} CA atoms. "
        f"Mode: {'shared fixed CA atoms' if align_on is True else 'explicit mapping (CA only)'}."
    )


class ContigJsonDataset(MolecularDataset):
    """
    Enables loading of JSON files containing contig data for benchmark design tasks,
    or the passing of examples through analogously-structured hydra configs.
    """

    def __init__(
        self,
        *,
        data: PathLike | Dict[str, dict | DesignInputSpecification],
        cif_parser_args: dict | None,
        transform: Transform | Compose | None,
        name: str | None,
        subset_to_keys: List[str] | None,
        eval_every_n: int,
    ):
        """
        Args:
            - data: path to the JSON file containing the contig data
            - cif_parser_args: arguments for the CIF parser
            - transform: transform to apply to the data
            - name: name of the dataset
            - subset_to_keys: list of keys to subset the data to
            - evaluate_every_n: how many times should this dataset be evaluated?
        """

        if isinstance(data, (PathLike, str)):
            self.json_path = data
            original_data = self._load_from_path(data)
        elif isinstance(data, DictConfig):
            self.json_path = None
            original_data = OmegaConf.to_object(data)
        else:
            self.json_path = None
            original_data = data

        # These will have already been added at inference time, but this block is useful for validation.
        if "global_args" in original_data:
            global_args = original_data.pop("global_args")
            for k, v in original_data.items():
                original_data[k].update(global_args)

        self._data = original_data

        if subset_to_keys is not None:
            assert (
                len(subset_to_keys) > 0
            ), "subset_to_keys must be a non-empty list of keys."
            self._data = {k: v for k, v in self._data.items() if k in subset_to_keys}
        self._check_json_keys()

        # ...basic assignments
        self.name = name if name is not None else "json-dataset"
        self.transform = transform

        self.cif_parser_args = cif_parser_args
        self.eval_every_n = eval_every_n

        if len(self) > 1_000:
            logger.warning(
                "ContigJsonDataset contains more than 1,000 entries. This may lead to performance issues."
            )
        elif len(self) == 0:
            raise ValueError(
                "ContigJsonDataset is empty, data: {}. Names: {}".format(
                    data, self.names
                )
            )

        l = 46
        fmt_names = textwrap.fill(
            ", ".join(self.names), width=l
        )  # .replace('\n', '+\n+ ')
        logger.info(
            f"\n+{l * '-'}+\n"
            f"Dataset {self.name}:\n"
            f"  - Found {len(self):,} examples:\n"
            f"{fmt_names}\n"
            f"\n+{l * '-'}+\n"
        )

    @staticmethod
    def _load_from_path(data):
        """Load data from a JSON or YAML file."""
        assert os.path.exists(data), f"Input file {data} does not exist."
        with open(data, "r") as f:
            if data.endswith(".json"):
                data = json.load(f)
            elif data.endswith(".yaml"):
                data = yaml.safe_load(f)
            else:
                raise ValueError(f"Input file {data} must be a JSON or YAML file.")
        return data

    def _check_json_keys(self):
        """Check if the JSON keys are valid."""
        for k, data in self.data.items():
            if not isinstance(data, (dict, DesignInputSpecification)):
                raise ValueError("Each item in the JSON data must be a dictionary.")

    @property
    def data(self):
        """Expose underlying dataframe as property to discourage changing it (can lead to unexpected behavior with torch ConcatDatasets)."""
        return self._data

    @property
    def names(self) -> List[str]:
        return list(self.data.keys())

    def __len__(self) -> int:
        """Pass through the length of the wrapped dataset."""
        return len(self.names)

    def __contains__(self, example_id: str) -> bool:
        """Pass through the contains method of the wrapped dataset."""
        return example_id in self.names

    def id_to_idx(self, example_id: str) -> int:
        """Pass through the id_to_idx method of the wrapped dataset."""
        return self.names.index(example_id)

    def idx_to_id(self, idx: int) -> str:
        """Pass through the idx_to_id method of the wrapped dataset."""
        return self.names[idx]

    def __getitem__(self, idx: int) -> Any:
        """Pass through the getitem method of the wrapped dataset."""
        example_id = self.idx_to_id(idx)
        spec = self.data[example_id]

        # if 'input' in metadata and not abspath, prepend the source json directory to the file path
        alt_spec_kwargs = None
        if not isinstance(spec, DesignInputSpecification):
            spec = ensure_input_is_abspath(spec, self.json_path)
            spec["cif_parser_args"] = self.cif_parser_args

            # Extract alt spec before creating the main DesignInputSpecification
            # (which has extra="forbid" and would reject unknown fields)
            alt_spec_kwargs = spec.get("alt", None)
            spec = {k: v for k, v in spec.items() if k != "alt"}

            spec = DesignInputSpecification.safe_init(**spec)

        # Create pipeline input
        data = spec.to_pipeline_input(example_id=example_id)

        # Build alt features if specified (separate spec, same pipeline)
        # Alignment must happen BEFORE transforms to preserve original residue numbering.
        if alt_spec_kwargs is not None:
            # Extract align_on before safe_init (DesignInputSpecification has extra="forbid")
            align_on = alt_spec_kwargs.pop("align_on", None)
            alt_spec_kwargs = ensure_input_is_abspath(alt_spec_kwargs, self.json_path)
            alt_spec_kwargs["cif_parser_args"] = self.cif_parser_args

            ## this block ensures that the number of designed tokens in the alt pipeline is the same as the main structure
            ## any number of designed tokens specified by the user in the alt contig input will be ignored,
            ## and substituted with the number on the main pipeline *at the end of the alt contig*
            ## NOTE: this block is designed in this way because it runs for every batch, so it needs to be robust even for the correct alt_spec_kwargs['contig']

            # strip alt contig of any designed portions
            stripped_contig = []
            for part in alt_spec_kwargs['contig'].split(','):
                if part == '/0' or any(c.isalpha() for c in part):
                    stripped_contig.append(part)
            alt_spec_kwargs['contig'] = ','.join(stripped_contig).strip(',/0')

            # add designed porstion equal to main contig
            for part in spec.extra['sampled_contig'].split(','):
                if not any(c.isalpha() for c in part) and part != '/0':
                    alt_spec_kwargs['contig'] += f',/0,{part}'
            
            ## end of block

            ## set as center of mass the same as what was picked for the main design
            com = data["specification"]["extra"]["com"]
            if len(com) != 3 or not isinstance(com[0], float):
                raise ValueError(f"Expected com to have only 3 float coordinates, got {com}")
            alt_spec_kwargs["ori_token"] = com


            alt_spec = DesignInputSpecification.safe_init(**alt_spec_kwargs)
            data_alt = alt_spec.to_pipeline_input(example_id=f"{example_id}_alt")

            # Align alt atom array to main BEFORE transforms (preserves original residue numbering)
            if align_on is not None:
                _align_alt_to_main(data, data_alt, align_on)

        # Apply transforms
        data = self.transform(data)

        if alt_spec_kwargs is not None:
            data_alt = self.transform(data_alt)
            # data["feats_alt"] = data_alt["feats"]
            # data["coord_atom_lvl_to_be_noised_alt"] = data_alt["coord_atom_lvl_to_be_noised"]
            data["alt"] = data_alt

        return data


def assemble_distributed_inference_loader_from_json(
    *, rank: int, world_size: int, **dataset_kwargs
) -> DataLoader:
    """
    Assemble a distributed inference DataLoader from JSONs.
    example:
        data={
            "backbone_0": {**args},
            "backbone_1": {**args}
        }
    """
    dataset = ContigJsonDataset(**dataset_kwargs)
    sampler = SequentialSampler(dataset)
    return assemble_distributed_loader(
        dataset=dataset,
        sampler=sampler,
        rank=rank,
        world_size=world_size,
    )
