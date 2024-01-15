

# asteRIa: Modeling Robust Interactions <img alt="asteRIa-logo-github2" src="https://github.com/marastadler/RobustInteractions/assets/61226497/6cdb3096-a3a4-4fe1-84be-21026cd7677b" align="right" width="200" />







This repository contains the code for the manuscript 'A statistical workflow for the robust detection of interaction effects between
chromatin modifications'.

We assume to be in a setting with an (binary) experimental design matrix L with n rows (observations) and q features. 
Under this design a corresponding multiple readout P with n rows and p columns is available.
In our case L contains q different chromatin modifications and P contains the binding strength of p Proteins to the n nucleosomes.
We are interested in the combinatorial (bi-order interaction) effect between those modifications on the binding behaviour of the proteins.




<img width="1602" alt="sketch_forwardmodel" src="https://user-images.githubusercontent.com/61226497/194323015-9c3c64cd-18a7-4c2e-80be-301b72998abe.png">




Our workflow includes the following main steps:

- a Lasso for hierarchical interactions with a strong hierarchy (Bien et al., 2013) to reduce number of spurious interaction effects.

- (Complementary pairs) stability selection for robust model selection in the data-scarce regime.

- Replicate consistency check to remove noisy data that might cause spurious interaction effects.


Moreover, we provide a script (`script_analysis_workflow.md`) that shows how downstream analyses and visalization of the results from our workflow can be performed. In this script all figures in the manuscript can be reproduced.


