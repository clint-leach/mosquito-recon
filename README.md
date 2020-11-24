# Linking mosquito surveillance to dengue fever through Bayesian mechanistic modeling

Leach CB, Hoeting JA, Pepin KM, Eiras AE, Hooten MB, Webb CT (2020) Linking mosquito surveillance to dengue fever through Bayesian mechanistic modeling. PLoS Negl Trop Dis 14(11): e0008868. https://doi.org/10.1371/journal.pntd.0008868

## Abstract

Our ability to effectively prevent the transmission of the dengue virus through targeted control of its vector, \emph{Aedes aegypti}, depends critically on our understanding of the link between mosquito abundance and human disease risk.
Mosquito and clinical surveillance data are widely collected, but linking them requires a modeling framework that accounts for the complex non-linear mechanisms involved in transmission.
Most critical are the bottleneck in transmission imposed by mosquito lifespan relative to the virus' extrinsic incubation period, and the dynamics of human immunity.
We developed a differential equation model of dengue transmission and embedded it in a Bayesian hierarchical framework that allowed us to estimate latent time series of mosquito demographic rates from mosquito trap counts and dengue case reports from the city of Vit\'oria, Brazil.
We used the fitted model to explore how the timing of a pulse of adult mosquito control influences its effect on the human disease burden in the following year.
We found that control was generally more effective when implemented in periods of relatively low mosquito mortality (when mosquito abundance was also generally low).
In particular, control implemented in early September (week 34 of the year) produced the largest reduction in predicted human case reports over the following year.
This highlights the potential long-term utility of broad, off-peak-season mosquito control in addition to existing, locally targeted within-season efforts.
Further, uncertainty in the effectiveness of control interventions was driven largely by posterior variation in the average mosquito mortality rate (closely tied to total mosquito abundance) with lower mosquito mortality generating systems more vulnerable to control.
Broadly, these correlations suggest that mosquito control is most effective in situations in which transmission is already limited by mosquito abundance.


## Repository

This repository holds the code, HMC and control simulation outputs, and manuscript draft.

The directories are:

* `Analysis`: scripts for processing data, running Stan, simulating control, and generating all the figures in the manuscript.
* `Code`: `.stan` source files containing the model specifications
* `Manuscript`: tex source files for PLOS NTD manuscript

## Data and HMC output

Data and model output can be downloaded from Figshare at https://doi.org/10.6084/m9.figshare.7905254

The code expects to find the data files (`.csv`) in a `Data` directory and model output (`.rds`) in a `Results` directory.
