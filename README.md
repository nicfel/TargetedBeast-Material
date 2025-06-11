# TargetedBeast-Material

This repository contains the data, analysis scripts, and figures associated with the publication:

**"Improving the Scalability of Bayesian Phylodynamic Inference through Efficient MCMC Proposals"**  
Remco R. Bouckaert, Paula H. Weidemüller, Luis R. Esquivel Gomez, and Nicola F. Müller

## Overview

Bayesian phylodynamic inference enables reconstruction of pathogen transmission history and evolutionary parameters using genome sequences. However, standard MCMC approaches struggle to scale beyond a few thousand taxa due to inefficient proposal distributions.

The **TargetedBeast** package introduces a set of novel MCMC operators for BEAST2 that:
- Target tree regions with higher mutation density using a pseudo-parsimony score
- Scale node heights more effectively using interval-based operators
- Increase mixing and reduce burn-in, especially for large datasets (e.g., 10,000+ sequences)

This repository contains:
- Scripts and code to reproduce all analyses and figures in the paper
- Benchmarks on simulated and real-world viral datasets
- Operator performance comparisons
- Supplementary figures and data

## Citation

If you use the code or concepts from this repository, please cite:

> Bouckaert RR, Weidemüller PH, Esquivel Gomez LR, Müller NF. Improving the Scalability of Bayesian Phylodynamic Inference through Efficient MCMC Proposals. (2025). [Manuscript in prep/submitted].

## License

All code and materials in this repository are released under the **LGPL 2.1** license.
