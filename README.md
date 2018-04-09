# TimemaConsOfAsex

This is the repository for the collected scripts used in the study "Consequences of asexuality in natural populations: insight from stick insects" (Bast et al. 2018, MBE).

The frozen code at the time of publication can be found here: https://github.com/jensbast/TimemaConsOfAsex/releases/tag/v1.0

### Components

- [pNpS estimates](pNpS) - scripts for calculation of synonymous and non-synonymous SNPs
- [dNdS estimates](dNdS) - script that uses alignments and RaXML branch-lengths calculate dN dS with codeml in the PAML package and outputs a table with all metrics
- [augment_OG.ipynb](augment_OG.ipynb) - script (ipython notebook) to mine ortholog variants from the OMA ouput
- [polymorphisms_processing.R](polymorphisms_processing.R) - R script to process output from pNpS program and calculate pS, pN and pNpS
- [commands_orth_algn_raxml](commands_orth_algn_raxml) - one-liner collection for getting the shared orthologous groups for all ten species, aligning them and running RAxML and preparation to run the dN/dS estimates
