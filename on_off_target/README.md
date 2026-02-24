
`pipeline.sh` runs the whole design pipeline:
1) RFD3 on the main on-target
2) ProteinMPNN on the structural samples, with the on-target, to generate final binder designs
3) RF3 on the designs with the on-target + metrics
4) RF3 on the designs with the off-target + metrics


The space of backbone conformations sampled when allowing the peptide backbone to move around is pretty ludicrous. I am not quite sure it makes a lot of sense. I would avoid it for the time being. Perhaps the sampling schedule can be tuned to be just right, but it's beyond the scope of what I want to do here.

RF3, in single-sequence mode, without templating, doesn't even know how to fold the MHC... not even the full one


DONE:
- compute re-folding metric: how well do design and folded structure agree? align MHCs and compute RMSD of pLDDT of the binders
- compute ipTM between binder and pMHC. Although the other components of ipTM are likely the same across designs since we have to use *exact* templating of the MHC... so maybe the global ipTM is good enough
- compute peptide-binder PAE
- compute amount of binder-peptide contact and amount of binder-MHC contact
- add proteinmpnn to the design pipeline

