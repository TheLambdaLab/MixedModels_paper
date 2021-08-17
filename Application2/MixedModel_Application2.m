function MixedModel_Application2(ESDOFstripes, MDOFstripes, EDPCapacity, EDP_col, plotflag)
% MixedModel_Application2(ESDOFstripes, MDOFstripes,EDPCapacity, EDP_col, plotflag)
% ----------------------------------------------------------------------------------
% INPUT
% ESDOFstripes: stripes of the ESDOF model
% MDOFstripes:  stripes of the MDOF model
% EDPCapacity:  vector with the drift capacities to be considered for the
%               mixed model
% EDP_col:      EDP that indicates collapse 
% plotflag       = if 1 the figures are plot for the results (default = 0)
%
% OUTPUT
% nothing but many figures are saved
% -----------------------------------------------------------------------------------
%% check inputs 
if nargin<4; plotflag=0; end
if nargin<3; error('needs at least 4 input arguments'); end

FigureFontSize=22; LabelSize=18;

%  range for which the fragility curves will be computed
IMrange = 0:0.001:10;
% ------------------------ NON-COLLAPSE FIT - MDOF MODEL ------------------
% fit the power-law model: EDP=aIM^b in the MDOF stripe(s)
% if only one stripe is available, then b=1 is adopted 

if length(MDOFstripes.IM) == 1
	% only one stripe is provided
	% the EDP=aIM^b passes from the median value of the stripe results
	% and has b=1
	beta_MDOF = 1;
	% if any nan appears in the EDP results, the median & std will return nan
	% in this case we will compute the median value as the 50% quantile
	if sum(isnan(MDOFstripes.EDP{1,1})) > 0
		med_EDP=quantile(MDOFstripes.EDP{1,1},0.5);
		sigma_MDOF=(quantile(log(MDOFstripes.EDP{1,1}),0.84)-quantile(log(MDOFstripes.EDP{1,1}),0.16))/2.;
	else
		med_EDP = median(MDOFstripes.EDP{1,1});
		sigma_MDOF = std(log(MDOFstripes.EDP{1,1}));
	end
	% med_EDP = a*MDOFstripes.IM^1 -->  alpha = med_EDP/MDOFstripes.IM
	alpha_MDOF = med_EDP/MDOFstripes.IM;
	
	% find the fraction of collapses in this stripe
	idx_noncol_MDOF = find(or(MDOFstripes.EDP{1,1} < EDP_col, ~isnan(MDOFstripes.EDP{1,1})));
	Pstripe = (length(MDOFstripes.EDP{1,1})-length(idx_noncol_MDOF))/length(MDOFstripes.EDP{1,1});
else
	% in this case the power law fit is employed to estimate the beta,
	% alpha and sigma parameters
	print('only one stripe is supported for now, slightly modify the code as per Example 2 so that it works for more than one stripes on the MDOFstripes')
end


% do the same for the ESDOF data
% in the non collapse fit we omit the stripes with more than
% 16% collapses since they will bias the non-collapse fit

% create a list with the IMs and EDPs
% initialize ESDOF IMs and EDPs that will be used in the non collapse fit
% also calculate the ratio of collapses per stripe that will be used in the
% collapse fit
IM_ESDOF_nc = [];
EDP_ESDOF_nc = [];
CollapseFraction_ESDOF = zeros(1, length(ESDOFstripes.IM));
Nruns = zeros(1, length(ESDOFstripes.IM));
Ncollapses = zeros(1, length(ESDOFstripes.IM));
idx_noncol =  cell(1, length(ESDOFstripes.IM));
for i=1:length(ESDOFstripes.IM)
	% find number of runs per stripe
	Nruns(1,i) = length(ESDOFstripes.EDP{i, 1});
	% find indices of no collapse points
	idx_noncol{1,i} = find(or(ESDOFstripes.EDP{i, 1} < EDP_col, ~isnan(ESDOFstripes.EDP{i, 1})));
	% find number of collapses in each stripe
	Ncollapses(1,i) = Nruns(1,i) - length(idx_noncol{1,i});
	CollapseFraction_ESDOF(1, i) = Ncollapses(1,i)/Nruns(1,i);
	if CollapseFraction_ESDOF(1, i) < 0.16
		IM_ESDOF_nc = [IM_ESDOF_nc ; ones(length(ESDOFstripes.EDP{i, 1})-Ncollapses(1,i),1)*ESDOFstripes.IM(i,1)];
		EDP_ESDOF_nc = [EDP_ESDOF_nc ; ESDOFstripes.EDP{i, 1}(idx_noncol{1,i})];
	end
end

% fit the power law model in the non collapse data of the ESDOF
[regress_ESDOF,~,~,~,stats_ESDOF]=regress(log(EDP_ESDOF_nc),[ones(length(IM_ESDOF_nc),1),log(IM_ESDOF_nc)]);
sigma_ESDOF = sqrt(stats_ESDOF(end));
alpha_ESDOF = exp(regress_ESDOF(1,1));
beta_ESOF = regress_ESDOF(2,1);

% fit the mle to estimate the collapse parameters
[b_ESDOF,d_ESDOF]=glmfit([log(ESDOFstripes.IM)],[Ncollapses',Nruns'],'binomial','link','probit');

% convert probit coefficients to lognormal distribution parameters
thetaC_ESDOF = exp(-b_ESDOF(1)/b_ESDOF(2));
betaC_ESDOF = 1/b_ESDOF(2);



% determine the alpha_Mixed, beta_Mixed and sigma_Mixed
% if P[collapse|IM] <= 0.16 then alpha_Mixed = alpha_MDOF
% else alpha_Mixed = alpha_ESDOF
if Pstripe <= 0.16
	alpha_Mixed = alpha_MDOF;
else
	alpha_Mixed = alpha_ESDOF;
end
beta_Mixed = beta_ESOF;
sigma_Mixed = sigma_ESDOF ;

% determine the theta and beta of the mixed model
if and(Pstripe >= 0.2, Pstripe <= 0.8)
	% eq. 8 from mixed model paper
	thetaC_Mixed = logninv(0.50, log(MDOFstripes.IM)-norminv(Pstripe)*betaC_ESDOF, betaC_ESDOF);
else
	thetaC_Mixed = thetaC_ESDOF;
end

if or(and(Pstripe > 0., Pstripe < 0.2), and(Pstripe > 0.8, Pstripe < 1.0))
	betaC_Mixed = (log(MDOFstripes.IM)-log(thetaC_ESDOF))/norminv(Pstripe);
else
	betaC_Mixed = betaC_ESDOF;
end



% for each limit state compute the fragility curves

% MDOF model
ProbNC_MDOF = zeros(length(EDPCapacity), length(IMrange));
% ESDOF model
ProbNC_ESDOF = zeros(length(EDPCapacity), length(IMrange));
ProbLS_ESDOF = zeros(length(EDPCapacity), length(IMrange));
%  Mixed model
ProbNC_Mixed = zeros(length(EDPCapacity), length(IMrange));
ProbLS_Mixed = zeros(length(EDPCapacity), length(IMrange));
% compute the collapse fragility curve of the ESDOF model
ProbC_ESDOF(1,:) = normcdf((log(IMrange/thetaC_ESDOF))/betaC_ESDOF);
% compute the collapse fragility of the mixed model
ProbC_Mixed(1,:) = normcdf((log(IMrange/thetaC_Mixed))/betaC_Mixed);
	
for i_edpc = 1:length(EDPCapacity)
	% compute the fragility curve for the MDOF model
	IMc50_MDOF=(EDPCapacity(1,i_edpc)/alpha_MDOF)^(1/beta_MDOF);
	ProbNC_MDOF(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_MDOF))/(sigma_MDOF/beta_MDOF));
	
	
	% ------ ESDOF -----
	% compute the non-collapse fragility curve for the ESDOF model
	IMc50_ESDOF=(EDPCapacity(1,i_edpc)/alpha_ESDOF)^(1/beta_ESOF);
	ProbNC_ESDOF(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_ESDOF))/(sigma_ESDOF/beta_ESOF));
	% compute the probability of exceeding the limit state for the ESDOF
	% model by combining the non-collapse with the collapse probability as
	% per the total probability theorem
	ProbLS_ESDOF(i_edpc,:)=ProbNC_ESDOF(i_edpc,:).*(ones(1,length(IMrange))-ProbC_ESDOF(1,:))+1.0*ProbC_ESDOF(1,:);
	
	% ------ Mixed -----
	% compute the probability of non-collapse for the mixed model
	% if P[collapse|IM] <= 0.16 on the MDOF model then alpha_MDOF is
	% adopted else alpha_ESDOF
	IMc50_Mixed = (EDPCapacity(1,i_edpc)/alpha_Mixed)^(1/beta_Mixed);
	ProbNC_Mixed(i_edpc,:)=normcdf((log(IMrange)-log(IMc50_Mixed))/(sigma_Mixed/beta_Mixed));
	% compute the fragility of the mixed model
	ProbLS_Mixed(i_edpc,:)=ProbNC_Mixed(i_edpc,:).*(ones(1,length(IMrange))-ProbC_Mixed(1,:))+1.0*ProbC_Mixed(1,:);
end
if plotflag == 1
	% plot non-collapse stripe analysis results for the NSP/ESDOF and the IDA/MDOF
	% models
	figure (1)
	hold on; grid on; box on
	% plot stripe analysis results of the ESDOF model
	stripes_ESDOF1_EDP = [];
	stripes_ESDOF1_IM = [];
	stripes_ESDOF2_EDP = [];
	stripes_ESDOF2_IM = [];
	for i=1:length(ESDOFstripes.IM)	
		EDPs = ESDOFstripes.EDP{i, 1};
		if CollapseFraction_ESDOF(1, i) > 0.16
			stripes_ESDOF2_EDP = [stripes_ESDOF2_EDP; EDPs(idx_noncol{1,i})];
			stripes_ESDOF2_IM = [stripes_ESDOF2_IM; ESDOFstripes.IM(i,1)*ones(length(idx_noncol{1,i}),1)];
		else
			stripes_ESDOF1_EDP = [stripes_ESDOF1_EDP; EDPs(idx_noncol{1,i})];
			stripes_ESDOF1_IM = [stripes_ESDOF1_IM; ESDOFstripes.IM(i,1)*ones(length(idx_noncol{1,i}),1)];
		end
	end
	plot(stripes_ESDOF1_EDP ,stripes_ESDOF1_IM, 'o', 'markerfacecolor', [128, 128, 128]/255, 'markeredgecolor', 'k', 'markersize', 6.0)
	plot(stripes_ESDOF2_EDP ,stripes_ESDOF2_IM, 'o', 'markerfacecolor',  [212, 208, 200]/255, 'markeredgecolor', 'k', 'markersize', 6.0)
	% plot stripe of the MDOF model
	EDPs_MDOF = MDOFstripes.EDP{1,1};
	ylim([0, 2.0])
	plot( EDPs_MDOF(idx_noncol_MDOF),MDOFstripes.IM*ones(length(idx_noncol_MDOF)), 'v', 'markerfacecolor', 'm', 'markersize', 10.0, 'markeredgecolor', 'k')
	set(gca,'fontsize',LabelSize)
	xlabel('maximum interstory drift ratio, \theta_{max}', 'fontsize', FigureFontSize)
	ylabel('Sa(T_1,5%) (g)', 'fontsize', FigureFontSize)
	legend('NSP/ESDOF (P[C|IM] < 16%)', 'NSP/ESDOF (P[C|IM] \geq 16%)', 'stripe MDOF')
	
	% plot figures for all limit states
	for i_edpc = 1:length(EDPCapacity)
		figure (i_edpc+1)
		hold on
		plot(IMrange, ProbLS_Mixed(i_edpc,:), 'color', [0,0.45,0.74], 'linewidth', 3)
		plot(IMrange, ProbLS_ESDOF(i_edpc,:), 'color', [0.31, 0.31, 0.31], 'linewidth', 3)
		xlabel('Sa(T_1,5%) (g)', 'fontsize', FigureFontSize)
		set(gca,'fontsize',LabelSize)
		hold on; grid on; box on
		ylabel(['P[\theta_{max}>\theta_{max,C}=', num2str(round(EDPCapacity(1,i_edpc),2)),'|Sa(T_1,5%)'], 'fontsize', FigureFontSize)
		legend('mixed model', 'NSP/ESDOF 5-parameter fit')
	end
end

end

