function [ EMD, CM] = DistributionDistance( dist1, dist2 )
% This function computes the absolute difference of two distributions and the
% Cramer-Von Mises distance is a simplified way, i.e. I sum the differences
% of the distributions.
% 
% INPUT
% dist1 = 1x3 cell array with: name of the inverse distribution, median value,
%         standard deviation. 
%         Alternatively, one can provide the P and x values directly for 
%         this distribution in cells 1,2 and 1,3, respectively. In this case we
%         don't need the name of the distribution thus the first entry
%         should be empty
% dist2 = 1x3 cell array with: name of the inverse distribution, median value,
%         standard deviation. 
%         Alternatively, one can provide the P and x values directly for this 
%         distribution in cells 1,2 and 1,3, respectively. In this case we
%         don't need the name of the distribution thus the first entry
%         should be empty
% OUTPUT
% EMD   = eaths moving distance of the distributions, i.e. the sum of the
%         absolute difference between the two cdfs
% CM    = Cramer-Von Mises distance, i.e. the square root of the sum of
%         squares of the differences of the two distributions

% Created 2019/12/8 by AC



%---- Example:
% dist1={'logncdf;',log(1.5),0.3};
% dist2={'logncdf',log(2.5),0.4};
%  [ EMD, CM] = DistributionDistance( dist1, dist2 )
%-----

if or(~length(dist1)==3,~length(dist1)==3); 
	error('input parameters should have 3 entries'); 
end
% check that dist1 and dist2 {1,2} and {1,3} inputs are of the same length
if or(~(length(dist1{1,2})==length(dist1{1,3})),~(length(dist2{1,2})==length(dist2{1,3}))); 
	error('check input variables'); 
end


% set the number of points to be created
N=10000;
% set the maximum x value if input parameters of dist1 and dist2
if and(length(dist1{2})~=1,length(dist2{2})~=1)
	xmax=min(max(dist1{3}),max(dist2{3}));
else 
	xmax=100;
end
% create N equally spaced probability values between 0 and 1
x=linspace( xmax./(2*N),xmax-xmax./(2*N), N);

% process distribution's 1 inputs
if length(dist1{2})~=1
	% check that dist1{1,2} input is probability values
	if or(min(dist1{2})<0,max(dist1{2})>1); error('dist1{1,2} should be probability values between 0 and 1'); end
	p1=interp1(dist1{3},dist1{2},x);
else 
	% get the probabilities corresponding to the x values given the
	% distribution and its parameters 
	p1=feval(dist1{1},x,dist1{2},dist1{3});
end


% process distribution's 2 inputs
if length(dist1{2})~=1
	% check that dist2{1,2} input is probabilty values
	if or(min(dist2{2})<0,max(dist2{2})>1); error('dist1{1,2} should be probability values between 0 and 1'); end
	p2=interp1(dist2{3},dist2{2},x);
else 
	% get the probabilities corresponding to the x values given the
	% distribution and its parameters
	p2=feval(dist2{1},x,dist2{2},dist2{3});
end


% find the differences of the probabilities for each x value
dif=p1-p2;
% absolute difference
EMD=sum(abs(dif));

% CM
CM=sqrt(sum(dif.^2));

end

