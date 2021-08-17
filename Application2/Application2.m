function Application2(bskdir, EDPcapacity, EDP_col, plotflag)
% Application2(bskdir, EDPcapacity, EDP_col, plotflag)
% calling this function you can run Application 2 of the "Mixed probabilistic seismic
% demand models for fragility assessment", A. Chatzidaki & D. Vamvatsikos,
% BEE. It reads stripe analyses results of the two source models and 
% computes the mixed model
% -----------------------------------------------------------------------------------
% INPUT
% bskdir:      basik directory where the Example2_MixedModel.m file is stored
% EDPcapacity: (1xNlimit states) array with the EDP capacities for the
%              limit states for which the fragility curves will be computed
% EDP_col:     EDP that indicates collapse 
% plotflag:    1 to plot the figures
% -----------------------------------------------------------------------------------
% OUTPUTS
% nothing but saves many files
% -----------------------------------------------------------------------------------
% Created by AC

% EXAMPLE:
% EDPcapacity = [0.01,0.02,0.03]
% Application2('.../Application_2',EDPcapacity , 0.08,1)


% the stripe analysis results of the NSP/ESDOF and the IDA/MDOF models are
% stored in the bskdir/StripeData path
InFilesPath = strcat(bskdir, '/StripeData');

% load the stripe analysis results of the two models
MDOFstripes = load(fullfile(InFilesPath,'IDA_MDOF.mat'));
ESDOFstripes = load(fullfile(InFilesPath,'NSP_ESDOF.mat'));

% run the main code
MixedModel_Application2(ESDOFstripes, MDOFstripes, EDPcapacity, EDP_col, plotflag)

end


