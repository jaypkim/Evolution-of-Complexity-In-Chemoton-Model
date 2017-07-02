# Evolution of Complexity in a Chemoton Model

This project focused on exploring some of the prerequisite conditions necessary for the formations of life in a modified chemoton model.
To do so, I used an agent-based model, using molecules as agents, to simulate the membrane-building biochemical pathway of a chemoton,
and used a performance scoring mechanism to track the chemoton's ability to expand its membrane and stably replicate, 
given external selective pressures. Here, the key modification made to the chemoton model is creating a gene-centric chemoton,
where the metabolic and membrane-building pathway is supplemented by protein products created from the "gene". The "gene" here should
not be confused for the modern form of a gene, but rather as an early RNA precursor.

## ABSTRACT

Due to recent advents in computing power and software, agent-based models (ABM) are becoming a popular method to simulate biological 
phenomena. However, much remains to be seen about its application on cellular processes at a mesoscopic level due to a lack of experimental data. The difficulty of parameter estimation, particularly for multivariate ABM models, also makes modeling multiple biological phenomena within a single model difficult. This project attempts to model a proto-cell on a mesoscopic level, using molecules as agents. I employed a novel approach to isolate a biochemical pathway of interest within a simulation without neglecting the effects of side reactions and flux in the proto-cell. I observed that environmental stress on membrane formation pathways, rate of gene expression, flux, and side-reactions may have non-linear impacts on the development of the proto-cell. Furthermore, the model used in this project 
demonstrates the effects of paralogous genes, gain-of-function mutations, and mutations that affect an enzyme’s capability on a 
biochemical pathway. Here, paralogous genes and gain-of-function mutations allow recovery of the proto-cell from otherwise fatal 
environmental conditions, whereas the greater variance from a greater likelihood and strength of mutations often kill off a large 
degree of the proto-cells.

## Using the code

This model was created in NetLogo 5.3.1, which is not the latest version of NetLogo. NetLogo does not support backwards compatability in its later iterations, so it is necessary to download 5.3.1 rather than the 6.0.1 version. To open and implement the model, use the .nlogo file. The .nls files are supplementary modules that are necessary to the model's functionality and must also be downloaded. The Enzyme_addition.txt file provides instructions on how to add a new type of agent (molecule) to the simulation.
