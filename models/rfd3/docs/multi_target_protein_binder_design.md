# RFdiffusion3 — Protein binder design with multiple targets

**NOTE**: this is intended as an extension of `ppi_design_tutorial.md` and `protein_binder_design.md`. Please familiarize with those before reading this.

## Overview

This README is for handling two use cases with RFD3:
1. Designing a protein binder for two targets in a single shot
2. Designing a protein binder for a target, while avoiding binding an off-target

We refer to both of these cases as Alternative-Target Guidance (ATG). ATG is handled in a similar way as Classifier-Free Guidance (CFG). While CFG follows a score that is the weighted sum of the target-conditioned score and the unconditioned score, ATG follows the weighted sum of two distinct target-conditioned scores, using different targets: `(1 - alt_scale) * score_t(x | main target) + alt_scale * score_t(x | alternative target)`.

If `0 < alt_scale < 1`, the score will design a binder for both targets with relative weight given by `alt_scale`. If `alt_scale < 0`, the score will design a binder to the main target, while avoiding the off target; in this case, `alt_scale` needs to be tuned carefully, as too extreme of a value might result in non-physical binders. A theoretical writeup of this procedure is coming soon.

**On aligning the main and alt structures.** A crucial consideration is that the diffusion trajectory is *not invariant* to the *relative position* of the main and alt structures. We can think of the design process as e.g. sampling from the probability distrution of binders conditioned on two targets *in a specific relative position* in space. Users can either provide the input structures in the desired relative orientation, or use the alignment specifications provided in the code and described below.


## How to specify the alternative target

To include the presence of an additional target, specify it in the field "alt", like in the following example:

```json
{
    "magea3": {
        "dialect": 2,
        "infer_ori_strategy": "hotspots",
        "input": "./input_pdbs/5brz_pmhc.pdb",
        "contig": "A1-275,/0,B1-100,/0,C1-9,/0,75-85",
        "select_hotspots": {
            "C1": "OE1,CD,OE2",
            "C4": "CB,CG,CD",
            "C5": "CG1,CD1",
            "C7": "CE1",
            "C8": "CD1"
        },
        "select_fixed_atoms": "A1-275,B1-100,C1-9",
        "is_non_loopy": true,
        "alt": {
            "dialect": 2,
            "input": "./input_pdbs/5bs0_pmhc.pdb",
            "contig": "A1-275,/0,B1-100,/0,C1-9",
            "select_fixed_atoms": "A1-275,B1-100,C1-9",
            "is_non_loopy": true,
            "align_on": {"A1-275,B1-100": "A1-275,B1-100"}
        }
    }
}
```
The specifications of the alt target go through *almost* the same pipeline as the main target. The same conditionings can be applied as those that can be applied to the main target (even hotspots, though they are not present in the above example).

Notable differences are as follows:
1. The "contig" field should only specify which parts of the alt "input" structure should be used as conditioning. Providing chains to design like for the main spec will result in an error.
2. "align_on", if present, prompts the code to align the main (key) and alt (value) "input" structures along the specified residues' CA atoms. In the example above, residues "A1-275,B1-100" of the main structure are aligned with the "A1-275,B1-100" residues of alt structure. While chains are resnums match in this example, they do not need to be; the only requirement is that the number of residues match. For example, a valid specification would be {"A10-22": "B35-47"}, but {"A10-22": "B35-48"} would be invalid and an error will be thrown.

The guidance parameter `alt_scale` can be set as a command-line argument to `rfd3 design`, e.g. `inference_sampler.alt_scale=0.5`, just like for the analogous CFG parameter.

