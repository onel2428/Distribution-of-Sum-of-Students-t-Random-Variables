function [v, s] = TFitting_2RVs_AbsoluteMoment(v1,v2,s1,s2)

%% Inputs
% [v1 v2] -> degrees of freedom of each of the two input t-RVs
% [s1 s2] -> scale of each of the two input t-RVs

%% Outputs
% v       -> degree of freedom of the fitting distribution 
% s       -> scale of the fitting distribution

% Relevant Moments of the Linear Combination
M2 =  SecondMoment([v1 v2],[s1 s2]);        % Second moment
MA = AbsoluteMoment_2RVs([v1 v2],[s1 s2]);  % Absolute Moment

% Fitting according to eqs. (32) and (33) in [ref]
v = (pi*MA-2*sqrt(2*M2))./(sqrt(pi)*MA-sqrt(2*M2));
s = sqrt(((pi-2*sqrt(pi))*MA*M2)./(pi*MA-2*sqrt(2*M2)));


%%
%[ref] Fitting the Distribution of Linear Combinations of $t-$Variables
%with more than 2 Degrees of Freedom (submitted to Journal of Probabilities and Statistics)

end