# Linking mosquito surveillance to dengue fever through Bayesian mechanistic modeling

This work was done in collaboration with Jennifer Hoeting, Kim Pepin, Alvaro Eiras, Mevin Hooten, and Colleen Webb.

## Abstract

Our ability to effectively prevent the transmission of the dengue virus through targeted control of its vector, Aedes aegypti, depends critically on our understanding of the link between mosquito abundance and human disease risk.
Mosquito and clinical surveillance data are widely collected, but linking them requires a modeling framework that accounts for the complex non-linear mechanisms involved in transmission.
Most critical are the bottleneck in transmission imposed by mosquito lifespan relative to the virus' extrinsic incubation period, and the dynamics of human immunity.
We developed a differential equation model of dengue transmission and embedded it in a Bayesian hierarchical framework that allowed us to estimate latent time series of mosquito demographic rates from mosquito trap counts and dengue case reports from the city of Vitoria, Brazil.
We used the fitted model to explore how the timing of a pulse of adult mosquito control influences its effect on the human disease burden in the following year.
We found that control implemented roughly 30 to 40 weeks prior to the start of the year produced the largest reduction in predicted human case reports for that year.
This highlights the potential long-term utility of broad, off-peak-season mosquito control in addition to existing, locally targeted within-season efforts.
Human immunity played an important role in mediating the effect of control, with the 30 to 40 week lag generally aligning with the estimated local minimum in the size of the human susceptible pool.
Further, uncertainty in the effectiveness of control interventions was driven largely by posterior variation in the period of cross-immunity and the average mosquito mortality rate, with shorter periods of cross-immunity and higher mosquito mortality leading to more effective control.
Understanding the dynamics of human immunity is thus critical for understanding how the effects of mosquito control propagate through the dengue system.

## Repository

This repository holds the code, HMC and control simulation outputs, and manuscript draft.

The directories are:

* Analysis: scripts for processing data, running Stan, simulating control, and generating all the figures in the manuscript.
* Code: `.stan` source files containing the model specifications
* Results: `.rds` files containing HMC chains (`gamma_eip_dv0.rds`) and output from mosquito control simulations (`gamma_control.rds`)

## Data

Data to which the model is fit can be download from Figshare and placed in a `Data` directory.
