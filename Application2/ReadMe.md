# **Application 2: Mixing different structural analysis types**
Using this set of files you can reproduce the results of Application 2 found in:
Mixed probabilistic seismic demand models for fragility assessment, 
Akrivi Chatzidaki & Dimiitrios Vamvatsikos, BEE
## Files contained
- StripeData:               folder where the stripe analysis results of the ESDOF and the MDOF models are stored
  * IDA_MDOF.mat: stripe analysis results of the MDOF model (7-records)
  * IDA_MDOF_44recs.mat: stripe analysis results of the MDOF model (44-records)
  * NSP_ESDOF.mat: stripe analysis results of the ESDOF model
- MixedModel_Application2.m:  fuction called when  running the Application2.m
- Application2.m:  main function needed for running the analysis
## How to run
Call the Application2 function via Matlab.
