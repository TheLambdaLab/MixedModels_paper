function Application1(basicdir, fiber_stripes_fname, lumped_stripes_fname, DOP_fname,EDPCapacity, plotflag)
% Application1(basicdir, fiber_stripes_fname, lumped_stripes_fname, DOP_fname,EDPCapacity, plotflag)
% calling this function you can run Application 1 of the "Mixed probabilistic seismic
% demand models for fragility assessment", A. Chatzidaki & D. Vamvatsikos, BEE
% It loads the stripe analysis results of the lumped and the distributed
% plasticity models and generates the mixed model
% -----------------------------------------------------------------------------------
% INPUT
% basicdir:    basic directory where the input data and the functions are stored
% fiber_stripes_fname: name of the .mat file containing the stripe analysis results 
%              of the fiber plasticity model. The .mat file should contain:
%              IM = (1xNstripes) array with the IM levels of the stripe analysis
%              EDP = {1xNstripes} cell with each entry being an (1xNruns) cell with 
%                    the stripe analysis results of the given stripe
%              desc = short description of the results contained within the file
% lumped_stripes_fname: name of the .mat file containing the stripe analysis results of
%              the lumped plasticity model with the same format as the 
%              fiber_stripes_fname variable
% DOP_fname:   name of the .mat file that contains the DOP per each
%              model, following the format:
%              confFiber: (1xNim) array with the DOP on the fiber model per confIM
%              confLumped: (1xNim) array with the DOP on the lumped model per confIM
%              confIM: (1xNim) array with the IM levels for which the DOP on each 
%                      model is determined 
%              The DOP per IM level should sum up to 1 when both models are considered. 
% EDPCapacity: (1xNedpc) array with the EDP capacities of the limit states
%              that will be considered in the mixed model
% plotflag:    if 1 the figures will be displayed (default: 1)
% -----------------------------------------------------------------------------------
% OUTPUT
% none but saves many figures in the bskdir/Results path
% -----------------------------------------------------------------------------------

% Example of application
% Application1('.../Application_1', 'Fiber_model.mat', 'Lumped_model.mat', 'DOP.mat',[0.015, 0.020, 0.025, 0.03], 1)
 
if nargin<6; plotflag=1; end

if plotflag==1
	% path where the results will be stored
	ResultsPath=[basicdir,'/Results'];
	% create this folder if it does not already exist
	if ~exist(ResultsPath, 'dir'); mkdir(ResultsPath);  end
end

% load stripe data of the two models
FiberModel = load(fullfile(basicdir, fiber_stripes_fname));
LumpedModel = load(fullfile(basicdir, lumped_stripes_fname));
% load DOP
DOP = load(fullfile(basicdir, DOP_fname));

% run the MixedModel_function function to find the mixed model
MixedModel_Application1(FiberModel, LumpedModel, DOP, EDPCapacity, 0.08, ResultsPath, plotflag)

end