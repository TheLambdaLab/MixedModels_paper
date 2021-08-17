function MixedModel_Application1( FiberModel, LumpedModel, DOP, EDPCapacity, EDP_col, ResultsPath, plotflag)
% MixedModel_Application1( FiberModel, LumpedModel, DOP, EDPCapacity, EDP_col, ...
%                         plotflag)
% -----------------------------------------------------------------------------------
% INPUT
% FiberModel =  structure with the stripe analysis results for the fiber model 
% LumpedModel = structure with the stripe analysis results for the lumped
%               plasticity model
% EDP_col =     EDP limit that indicates collapse 
% DOP =         structure with the DOP per IM level and model
% EDPCapacity = vector with the drift capacities you want to examine.
%               The mixed model will be fitted considering all
%               capacities defined herein
% ResultsPath  = path where the results will be stored (only if plotflag=1)
% plotflag     = if 1 the figures are plot for the results (default = 0)
% -----------------------------------------------------------------------------------
% OUTPUT
% many figures are generated
% -----------------------------------------------------------------------------------

if nargin <7; plotflag=0; end


% range of the IM for which the fragilities are computed
IMrange = 0:0.001:10.0;
% check that the DOPs for each IM level are defined for all models
if ~and(length(DOP.confFiber)== length(DOP.confLumped),length(DOP.confLumped)== length(DOP.confIM));
	error('check that a DOP value is provided for each IM')
end
% check that the DOPs of all models at each IM level sum up to 1.0
if DOP.confFiber+DOP.confLumped ~= ones(1, length(DOP.confFiber))
	error('DOP of the two models should sum up to one for each IM level')
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% STEP 3: fit the 5-parameter surrogate to source model i 
% ------------------------ FiberModel: ------------------------------------
% in this in the non collapse fit we will omit the stripes with more than
% 16% collapses since they will bias the non-collapse fit

% create a list with the IMs and EDPs that will be used in the non-collapse
% fit and also calculate the ratio of collapses per stripe that will be used in the
% collapse fit
IM_FiberModel_nc = [];
EDP_FiberModel_nc = [];
CollapseFraction_FiberModel = zeros(1, length(FiberModel.IM));
Nruns = zeros(1, length(FiberModel.IM));
Ncollapses = zeros(1, length(FiberModel.IM));
idx_noncol =  cell(1, length(FiberModel.IM));
for i=1:length(FiberModel.IM)
	% find number of runs per stripe
	Nruns(1,i) = length(FiberModel.EDP{i, 1});
	% find indices of non collapse points
	idx_noncol{1,i} = find(or(FiberModel.EDP{i, 1} < EDP_col, ~isnan(FiberModel.EDP{i, 1})));
	% find number of collapses in each stripe
	Ncollapses(1,i) = Nruns(1,i) - length(idx_noncol{1,i});
	CollapseFraction_FiberModel(1, i) = Ncollapses(1,i)/Nruns(1,i);
	if CollapseFraction_FiberModel(1, i) < 0.16
		IM_FiberModel_nc = [IM_FiberModel_nc ; ones(length(FiberModel.EDP{i, 1})-Ncollapses(1,i),1)*FiberModel.IM(i,1)];
		EDP_FiberModel_nc = [EDP_FiberModel_nc ; FiberModel.EDP{i, 1}(idx_noncol{1,i})];
	end
end


if length(unique(IM_FiberModel_nc)) == 1
	% only one stripe is provided
	% the EDP=aIM^b passes from the median value of the stripe results
	% and has b=1
	beta_FiberModel = 1;
	% if any nan appears in the EDP results, the median & std will return nan
	% in this case we will compute median and std via quantiles
	if sum(isnan(EDP_FiberModel_nc)) > 0
		med_EDP=quantile(EDP_FiberModel_nc,0.5);
		sigma_MDOF=(quantile(log(EDP_FiberModel_nc),0.84)-quantile(log(EDP_FiberModel_nc),0.16))/2.;
	else
		med_EDP = median(EDP_FiberModel_nc);
		sigma_MDOF = std(log(EDP_FiberModel_nc));
	end
	% med_EDP = alpha*MDOFstripes.IM^1 -->  alpha = med_EDP/MDOFstripes.IM
	alpha_FiberModel = med_EDP/FiberModel.IM;
else
	if length(EDP_FiberModel_nc) ~= 0
		% fit the power law model in the no collapse data of the ESDOF
		[regress_FiberModel,~,~,~,stats_FiberModel]=regress(log(EDP_FiberModel_nc),[ones(length(IM_FiberModel_nc),1),log(IM_FiberModel_nc)]);
		sigma_FiberModel = sqrt(stats_FiberModel(end));
		alpha_FiberModel = exp(regress_FiberModel(1,1));
		beta_FiberModel = regress_FiberModel(2,1);
	else
		error('no non-collapse points found for the fiber model')
	end
end
% compute non-collapse fragilities for all limit states of interest
% initialize probability of non collapse
ProbNC_Fiber = zeros(length(EDPCapacity), length(IMrange));
for i_edpc = 1:length(EDPCapacity)
	IMc50_Fiber=(EDPCapacity(1,i_edpc)/alpha_FiberModel)^(1/beta_FiberModel);
	ProbNC_Fiber(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_Fiber))/(sigma_FiberModel/beta_FiberModel));
end



% collapse fit & collapse fragility
% initialize probability of collapse
ProbC_Fiber = zeros(length(EDPCapacity), length(IMrange));

if sum(CollapseFraction_FiberModel) ~= 0
	% fit the mle to estimate the collapse parameters
	[b_FiberModel,d_FiberModel]=glmfit([log(FiberModel.IM)],[Ncollapses',Nruns'],'binomial','link','probit');
	% convert probit coefficients to lognormal distribution parameters
	thetaC_FiberModel = exp(-b_FiberModel(1)/b_FiberModel(2));
	betaC_FiberModel = 1/b_FiberModel(2);
	% compute collapse fragility
	ProbC_Fiber(i_edpc,:) = normcdf((log(IMrange/thetaC_FiberModel))/betaC_FiberModel);
end

% fragility of limit state
% initialize fragilities
ProbLS_Fiber = zeros(length(EDPCapacity), length(IMrange));
for i_edpc = 1:length(EDPCapacity)
	ProbLS_Fiber(i_edpc,:)=ProbNC_Fiber(i_edpc,:).*(ones(1,length(IMrange))-ProbC_Fiber(i_edpc,:))+1.0*ProbC_Fiber(i_edpc,:);
end



% ------------------------ LumpedModel: ------------------------------------
% create a list with the IMs and EDPs that will be used in the non-collapse
% fit and also calculate the ratio of collapses per stripe that will be used in the
% collapse fit
IM_LumpedModel_nc = [];
EDP_LumpedModel_nc = [];
CollapseFraction_LumpedModel = zeros(1, length(LumpedModel.IM));
Nruns = zeros(1, length(LumpedModel.IM));
Ncollapses = zeros(1, length(LumpedModel.IM));
idx_noncol =  cell(1, length(LumpedModel.IM));
for i=1:length(LumpedModel.IM)
	% find number of runs per stripe
	Nruns(1,i) = length(LumpedModel.EDP{i, 1});
	% find indices of no collapse points
	idx_noncol{1,i} = find(or(LumpedModel.EDP{i, 1} < EDP_col, ~isnan(LumpedModel.EDP{i, 1})));
	% find number of collapses in each stripe
	Ncollapses(1,i) = Nruns(1,i) - length(idx_noncol{1,i});
	CollapseFraction_LumpedModel(1, i) = Ncollapses(1,i)/Nruns(1,i);
	if CollapseFraction_LumpedModel(1, i) < 0.16
		IM_LumpedModel_nc = [IM_LumpedModel_nc ; ones(length(LumpedModel.EDP{i, 1})-Ncollapses(1,i),1)*LumpedModel.IM(i,1)];
		EDP_LumpedModel_nc = [EDP_LumpedModel_nc ; LumpedModel.EDP{i, 1}(idx_noncol{1,i})];
	end
end


if length(unique(IM_LumpedModel_nc)) == 1
	% only one stripe is provided
	% the EDP=aIM^b passes from the median value of the stripe results
	% and has b=1
	beta_LumpedModel = 1;
	% if any nan appears in the EDP results, the median & std will return nan
	% in this case we will compute the median value as the 50% quantile
	if sum(isnan(EDP_LumpedModel_nc)) > 0
		med_EDP=quantile(EDP_LumpedModel_nc,0.5);
		sigma_LumpedModel=(quantile(log(EDP_LumpedModel_nc),0.84)-quantile(log(EDP_LumpedModel_nc),0.16))/2.;
	else
		med_EDP = median(EDP_LumpedModel_nc);
		sigma_LumpedModel = std(log(EDP_LumpedModel_nc));
	end
	% med_EDP = alpha*MDOFstripes.IM^1 -->  alpha = med_EDP/MDOFstripes.IM
	alpha_LumpedModel = med_EDP/unique(IM_LumpedModel_nc);
else
	if length(EDP_LumpedModel_nc) ~= 0
		% fit the power law model in the no collapse data of the ESDOF
		[regress_LumpedModel,~,~,~,stats_LumpedModel]=regress(log(EDP_LumpedModel_nc),[ones(length(IM_LumpedModel_nc),1),log(IM_LumpedModel_nc)]);
		sigma_LumpedModel = sqrt(stats_LumpedModel(end));
		alpha_LumpedModel = exp(regress_LumpedModel(1,1));
		beta_LumpedModel = regress_LumpedModel(2,1);
	else
		error('no non-collapse points found for the fiber model')
	end
end
% compute non-collapse fragilities for all limit states of interest
% initialize probability of non-collapse
ProbNC_Lumped = zeros(length(EDPCapacity), length(IMrange));
for i_edpc = 1:length(EDPCapacity)
	IMc50_Lumped=(EDPCapacity(1,i_edpc)/alpha_LumpedModel)^(1/beta_LumpedModel);
	ProbNC_Lumped(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_Lumped))/(sigma_LumpedModel/beta_LumpedModel));
end

% collapse fit
%  initialize probability of collapse
ProbC_Lumped = zeros(1, length(IMrange));
if sum(CollapseFraction_LumpedModel) == 0
	% not any collapse point is available in the stripe thus the
	% probability of collapse is assumed to be zero
	pass
else
	% fit the mle to estimate the collapse parameters
	[b_LumpedModel,~]=glmfit([log(LumpedModel.IM)],[Ncollapses',Nruns'],'binomial','link','probit');
	% convert probit coefficients to lognormal distribution parameters
	thetaC_LumpedModel = exp(-b_LumpedModel(1)/b_LumpedModel(2));
	betaC_LumpedModel = 1/b_LumpedModel(2);
	% compute collapse fragility
	ProbC_Lumped = normcdf((log(IMrange/thetaC_LumpedModel))/betaC_LumpedModel);
end

% fragility of limit state
% initialize fragilities
ProbLS_Lumped = zeros(length(EDPCapacity), length(IMrange));
for i_edpc = 1:length(EDPCapacity)
	ProbLS_Lumped(i_edpc,:)=ProbNC_Lumped(i_edpc,:).*(ones(1,length(IMrange))-ProbC_Lumped)+1.0*ProbC_Lumped;
end




% pre-process DOP
% modify the confIM to add 0 in the first entry if it does not exist
% and a quite large value for the maximum value (if DOP is defined in ascending 
% order).
if sum(diff(DOP.confIM)>0)==length(DOP.confIM)-1; % this means that the DOP is set in ascending order
	if DOP.confIM(1)~=0;
		DOP.confIM=[0, DOP.confIM];
		DOP.confFiber = [DOP.confFiber(1), DOP.confFiber];
		DOP.confLumped = [DOP.confLumped(1), DOP.confLumped];
	elseif DOP.confIM(end)~=max(IMrange);
		DOP.confIM=[DOP.confIM, max(IMrange)];
		DOP.confFiber = [DOP.confFiber, DOP.confFiber(end)];
		DOP.confLumped = [DOP.confLumped, DOP.confLumped(end)];
	end
elseif sum(diff(DOP.confIM)<0)==length(DOP.confIM)-1; % this means that the DOP is set in descending order
	if DOP.confIM(1)~=max(IMrange);
		DOP.confIM=[max(IMrange), DOP.confIM];
		DOP.confFiber = [DOP.confFiber(1), DOP.confFiber];
		DOP.confLumped = [DOP.confLumped(1), DOP.confLumped];
	elseif DOP.confIM(end)~=0;
		DOP.confIM=[DOP.confIM,0];
		DOP.confFiber = [DOP.confFiber, DOP.confFiber(end)];
		DOP.confLumped = [DOP.confLumped, DOP.confLumped(end)];
	end
else
	error('confIM should be monotonically increasing or decreasing')
end


% interpolate the DOP of each model to get a DOP value for
% each IM in the IMrange
% initialize cells
confFiber=interp1(DOP.confIM,DOP.confFiber,IMrange);
confLumped=interp1(DOP.confIM,DOP.confLumped,IMrange);

% compute the target fragilities for each limit state by combining the
% fragilities of each limit state
% initialize target fragility
ProbLS_Tgt=zeros(length(EDPCapacity),length(IMrange));
for i_edpc=1:length(EDPCapacity)
	ProbLS_Tgt(i_edpc,:) = confFiber.*ProbLS_Fiber(i_edpc,:)+confLumped.*ProbLS_Lumped(i_edpc,:);
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% iterate on weight and compute fragilities 
WeightRange = 0:1:100;
Nweights=length(WeightRange);


% intialize cells
FiberWeight=zeros(1,Nweights);
LumpedWeight=zeros(1,Nweights);

weightNC_Lumped=cell(1,Nweights);
weightNC_Fiber = cell(1,Nweights);

ProbNC_Mixed=cell(1,Nweights);
ProbLS_Mixed=cell(1,Nweights);
EMDtgt_mixed=zeros(length(EDPCapacity),Nweights);
CMtgt_mixed=zeros(length(EDPCapacity),Nweights);
EMDfiber_mixed=zeros(length(EDPCapacity),Nweights);
EMDmixed_lumped=zeros(length(EDPCapacity),Nweights);

regress_Mixed=cell(1,Nweights);
stats_Mixed=cell(1,Nweights);
alpha_Mixed=cell(1,Nweights);
beta_Mixed=cell(1,Nweights);
sigma_Mixed=cell(1,Nweights);
for kk=1:length(WeightRange)
	% find the weight on each model. The relative weights should sum
	% up to one
	FiberWeight(1,kk) = WeightRange(kk)/100;
	LumpedWeight(1,kk) = 1-FiberWeight(1,kk); 

	% ~~~~~~~~~~~~~~~~~~~~~~~~ COLLAPSE FIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	% the collapse fit is directly coming from the lumped model
	
	% ~~~~~~~~~~~~~~~~~~~~~~~~ NON-COLLAPSE FIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	% employ weighted linear regression to determine the power-law fit of
	% the non-collapse data
	
	% i divide the weights by two on each stripe.
	if LumpedWeight(1,kk)==1.0
		% the non-collapse fit is determined from the lumped model
		sigma_Mixed{1,kk} = sigma_LumpedModel;
		alpha_Mixed{1,kk} = alpha_LumpedModel;
		beta_Mixed{1,kk} = beta_LumpedModel;
	elseif FiberWeight(1,kk)==1.0
		% the non-collapse fit is determined from the fiber model
		sigma_Mixed{1,kk} = sigma_FiberModel;
		alpha_Mixed{1,kk} = alpha_FiberModel;
		beta_Mixed{1,kk} = beta_FiberModel;
	else
		% weighted linear regression
		% find how many non collapse points are available for both models
		NNC_points=length(EDP_LumpedModel_nc)+length(EDP_FiberModel_nc);
		
		% compute the weights of the non-collapse fiber data
		weightNC_Fiber{1,kk} = NNC_points/length(EDP_FiberModel_nc)*FiberWeight(1,kk);
		weightNC_Lumped{1,kk} = NNC_points/length(EDP_LumpedModel_nc)*LumpedWeight(1,kk);
		weightsNC=[ones(length(EDP_FiberModel_nc),1)*weightNC_Fiber{1,kk};ones(length(EDP_LumpedModel_nc),1)*weightNC_Lumped{1,kk}];
		NCims=[IM_FiberModel_nc; IM_LumpedModel_nc];
		NCedps=[EDP_FiberModel_nc; EDP_LumpedModel_nc];
		if ~sum(weightsNC)==NNC_points
			warning('check the weights assigned on each model!')
		end
		[regress_Mixed{1,kk},~,~,~,stats_Mixed{1,kk}]=nlinfit(log(NCims),log(NCedps),@powerlaw,[1,1],'weights',weightsNC);
		sigma_Mixed{1,kk} = sqrt(stats_Mixed{1,kk}(end));
		alpha_Mixed{1,kk} = exp(regress_Mixed{1,kk}(1,1));
		beta_Mixed{1,kk} = regress_Mixed{1,kk}(1,2);
	end
	% calculate the fragility curve of the mixed model
	
	% mixed model
	for i_edpc = 1:length(EDPCapacity)
		% non-collapse data
		IMc50_Mixed=(EDPCapacity(1,i_edpc)/alpha_Mixed{1,kk})^(1/beta_Mixed{1,kk});
		ProbNC_Mixed{1,kk}(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_Mixed))/(sigma_Mixed{1,kk}/beta_Mixed{1,kk}));
		
		% collapse data coming from the lumped model
		ProbLS_Mixed{1,kk}(i_edpc,:)=ProbNC_Mixed{1,kk}(i_edpc,:).*(ones(1,length(IMrange))-ProbC_Lumped)+1.0*ProbC_Lumped;
		
		% compute the distance of the mixed model fragility from the
		% target fragility of the given limit state
		[EMDtgt_mixed(i_edpc,kk),CMtgt_mixed(i_edpc,kk)]=DistributionDistance( {[],ProbLS_Tgt(i_edpc,:)',IMrange}, {[],ProbLS_Mixed{1,kk}(i_edpc,:),IMrange});
		% find the distance of the distribution from the parent distributions
		[EMDfiber_mixed(i_edpc,kk),~]=DistributionDistance( {[],ProbLS_Fiber(i_edpc,:)',IMrange}, {[],ProbLS_Mixed{1,kk}(i_edpc,:)',IMrange});
		[EMDmixed_lumped(i_edpc,kk),~]=DistributionDistance( {[],ProbLS_Lumped(i_edpc,:)',IMrange}, {[],ProbLS_Mixed{1,kk}(i_edpc,:)',IMrange});
		% caclulate the distance between the lumped and the fiber model
		[EMDlumped_lumped,~]=DistributionDistance( {[],ProbLS_Lumped(i_edpc,:)',IMrange}, {[],ProbLS_Fiber(i_edpc,:)',IMrange});
		% if for any mm,kk pair the distance of mod1-mixed or mod2-mixed
		% is greater than this of mod1-mod2 then the EMD, CM and K
		% distances are "penalized" in order not to select them so I add an
		% extremely high value of e.g. 1000000000000
		if or(EMDlumped_lumped<EMDfiber_mixed(i_edpc,kk),EMDlumped_lumped<EMDmixed_lumped(i_edpc,kk))
			EMDtgt_mixed(i_edpc,kk)=EMDtgt_mixed(i_edpc,kk)+1000000000000;
		end
		
	end % limit state loop
	
	
	
end; % kk=1:length(WeightRange)



display('finding the weights that minimize the cumulative distance')

% find the weight that minimizes the sum of the distances for all the drift
% capacities examined
[~,indEMD]=min(sum(EMDtgt_mixed));


% now plot the results for this weight
if plotflag==1
	figure
	plot(DOP.confFiber,DOP.confIM,'k-','linewidth',3)
	hold on
	plot(DOP.confLumped,DOP.confIM,'k--','linewidth',3)
	xlabel('DOP','fontsize',22)
	ylabel('AvgSa (g)','fontsize',22)
	legend('fiber','lumped')
	set(gca,'fontsize',18)
	ylim([0, 1.5])
	grid on; box on
	saveas(gca,fullfile(ResultsPath,'DOP.fig'));
	saveas(gca,fullfile(ResultsPath,'DOP.epsc'), 'epsc2');
	saveas(gca,fullfile(ResultsPath,'DOP.emf'));
	saveas(gca,fullfile(ResultsPath,'DOP.png'));
	for i_edpc=1:length(EDPCapacity)
		figure
		hold on
		plot(IMrange,ProbLS_Fiber(i_edpc,:)','k-','linewidth',3)
		plot(IMrange,ProbLS_Lumped(i_edpc,:)','k--','linewidth',3)
		plot(IMrange,ProbLS_Tgt(i_edpc,:)','color',[128/254, 128/254, 128/254],'linewidth',3)
		plot(IMrange,ProbLS_Mixed{1,indEMD}(i_edpc,:)','m','linewidth',3)
		
		display(strcat('fiber, edpc= ', num2str(EDPCapacity(1,i_edpc))))
		vals = interp1(ProbLS_Fiber(i_edpc,50:500),IMrange(50:500),[0.16,0.50,0.84]);
		display(strcat('median= ', num2str(vals(2))))
		display(strcat('std= ', num2str(0.5*(log(vals(end))-log(vals(1))))))
		display(strcat('lumped, edpc= ', num2str(EDPCapacity(1,i_edpc))))
		vals = interp1(ProbLS_Lumped(i_edpc,50:1500),IMrange(50:1500),[0.16,0.50,0.84]);
		display(strcat('median= ', num2str(vals(2))))
		display(strcat('std= ', num2str(0.5*(log(vals(end))-log(vals(1))))))
		display(strcat('mixed, edpc= ', num2str(EDPCapacity(1,i_edpc))))
		vals = interp1(ProbLS_Mixed{1,indEMD}(i_edpc,50:1500),IMrange(50:1500),[0.16,0.50,0.84]);
		display(strcat('median= ', num2str(vals(2))))
		display(strcat('std= ', num2str(0.5*(log(vals(end))-log(vals(1))))))
		
		
		
		grid on; box on;
		set(gca,'fontsize',18)
		xlim([0 1.5])
		xlabel('AvgSa (g)','fontsize',22)
		ylabel(['P[{\theta}_{max}>',num2str(EDPCapacity(i_edpc)*100),'%|AvgSa]'],'fontsize',22)
		legend('fiber','lumped','target','mixed')
		saveas(gca,fullfile(ResultsPath,['Fragilities_LS',num2str(i_edpc),'.fig']));
		saveas(gca,fullfile(ResultsPath,['Fragilities_LS',num2str(i_edpc),'.epsc']), 'epsc2');
		saveas(gca,fullfile(ResultsPath,['Fragilities_LS',num2str(i_edpc),'.emf']));
		saveas(gca,fullfile(ResultsPath,['Fragilities_LS',num2str(i_edpc),'.png']));
	end
	
	
	% plot the non-collapse stripe analysis results of both models
	figure
	plot(EDP_LumpedModel_nc, IM_LumpedModel_nc, 'd', 'markersize', 10, 'markeredgecolor', 'k', 'markerfacecolor', [204/254, 204/254, 204/254], 'linewidth', 2)
	hold on
	plot(EDP_FiberModel_nc, IM_FiberModel_nc, 'o', 'markersize', 10, 'markeredgecolor', 'k', 'linewidth', 2)
	legend('lumped (NC)','fiber (NC)')
	grid on; box on;
	set(gca,'fontsize',18)
	xlim([0 0.07])
	ylim([0 1.5])
	xlabel('maximum interstory drift ratio, \theta_{max}','fontsize',22)
	ylabel('AvgSa (g)','fontsize',22)
	saveas(gca,fullfile(ResultsPath,'NC_Stripe_Data.fig'));
	saveas(gca,fullfile(ResultsPath,'NC_Stripe_Data.epsc'), 'epsc2');
	saveas(gca,fullfile(ResultsPath,'NC_Stripe_Data.emf'));
	saveas(gca,fullfile(ResultsPath,'NC_Stripe_Data.png'));
	% do the same but also add the fits
	figure
	loglog(EDP_LumpedModel_nc, IM_LumpedModel_nc, 'd', 'markersize', 10, 'markeredgecolor', 'k', 'markerfacecolor', [204/254, 204/254, 204/254], 'linewidth', 2)
	hold on
	plot(EDP_FiberModel_nc, IM_FiberModel_nc, 'o', 'markersize', 10, 'markeredgecolor', 'k', 'linewidth', 2)
	FiberEDP=alpha_FiberModel*IMrange.^beta_FiberModel;
	LumpedEDP=alpha_LumpedModel*IMrange.^beta_LumpedModel;
	MixedEDP=alpha_Mixed{1,indEMD}*IMrange.^beta_Mixed{1,indEMD};
	plot(FiberEDP, IMrange, 'k-', 'linewidth', 3)
	plot(LumpedEDP, IMrange, 'k--', 'linewidth', 3)
	plot(MixedEDP, IMrange, 'm', 'linewidth', 3)
	xlim([10^-4, 10^0])
	ylim([10^-2, 1.5])
	legend('lumped (NC)','fiber (NC)', 'fiber power-law fit','lumped power-law fit','mixed power law fit' )
	grid on; box on;
	set(gca,'fontsize',18)
	xlabel('maximum interstory drift ratio, \theta_{max}','fontsize',22)
	ylabel('AvgSa (g)','fontsize',22)
	saveas(gca,fullfile(ResultsPath,'NC_fit.fig'));
	saveas(gca,fullfile(ResultsPath,'NC_fit.epsc'), 'epsc2');
	saveas(gca,fullfile(ResultsPath,'NC_fit.emf'));
	saveas(gca,fullfile(ResultsPath,'NC_fit.png'));
end



end

