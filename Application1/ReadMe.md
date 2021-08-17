# Application 1: Mixing structural models of different accuracy
Using this set of files you can reproduce the results of Application 1 found in: Mixed probabilistic seismic demand models for fragility assessment, Akrivi Chatzidaki & Dimiitrios Vamvatsikos, BEE

## Files contained
- Lumped_model.mat: stripe analysis results of the lumped model (9 records per stripe/2 stripes)
- Fiber_model.mat: stripe analysis results of the fiber model (9 records per stripe/2 stripes)
- DOP.mat: degree of preference on each model
- MixedModel_Application1.m: function called when running the Application1.m
- Application1.m: main function needed for running the analysis
- powerlaw.m: function called for fitting the power law model
- DistributionDistance.m: function needed for computing the distance of two distributions
## How to run
Call the Application1 function via Matlab.
