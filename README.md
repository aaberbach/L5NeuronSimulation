# L5NeuronSimulation

Analyzing the firing properties of a single, detailed L5b PC neuron.

## Folder Structure

### Temporary
## These folders might be irrelevant or should be removed before publishing.

Attenuation - investigated the attenuation of the cell.

Base - what used to be the core simulation. PROBABLY REMOVE

EfficientBuild - merely a folder I used to experiment with making the build step use less memory. PROBABLY REMOVE

ErrorSynapse - used to investigate previous pyr2pyr.mod file. FIXED, REMOVE

F-ICurve - used to plot F-ICurve. NOT UP TO DATE

Gamma-rhythmic - gamma-rhythmic from a long time ago. OUTDATED

Group_EPSP_Tuning - excitatory tuning from when we were using Song et al. EPSP distribution. UNECESSARY

IPSP_Tuning - was going to be used for tuning Apical inhibiton, but we decided to just use IPSCs from basal inhibition to tune instead. UNECESSARY

Old_Folders - where I put some outdated folders like the ones above. PROBABLY REMOVE

Resistance - looking at input resistance.

### Core Folders

biophys_components - this is where mechanisms and templates are stored. This folder is used as the basis for all simulations in the project.

Clustering - this is where I worked on building the clustering mechanisms. MAYBE REMOVE

EPSP_Tuning - this folder has both old EPSP tuning for Song et al. and also the conductance scaling MAYBE REFACTOR

FullSimulation - this is the main folder, where the full simulations are run. See README in folder for more.

Group_EPSC_Tuning - where EPSC tuning is performed (on both basal and apical dendrites).

HayFiringProperties - confirming that our cell model matches the experimental firing properties that Hay et al. 2011 tuned for.

IPSC_Tuning - where IPSC tuning is performed (perisomatic, basal, and apical).

MorphAnalysis - folder for exploring the morphology of the cell.

Plot_Cell - folder for plotting the cell and synapse locations on the cell.

ReleaseProbability - used to build and check the release probability. COULD REMOVE

ShortPlasticity - used to build and check short term plasticity for both excitation and inhibition.

StepCurrent - used to ensure that a suprathreshold step current induces the proper behavior (as determined by Hay et al. 2011).

