function [v, s] = TFitting_KRVs_AbsoluteMoment(v,s)

%% Inputs
% v -> vector with the degrees of freedom of each of the input t-RVs
% s -> vector with the scale of each of the input t-RVs

%% Outputs
% v  -> degree of freedom of the fitting distribution 
% s  -> scale of the fitting distribution

%% Iterative procedure that repeatedly uses the fitting of the sum of two t-RVs as described in Section 3.1.1 of [ref]
v1 = v(1);
v2 = v(2);
s1 = s(1);
s2 = s(2);

for i=3:length(v)
    [v1, s1] = TFitting_2RVs(v1,v2,s1,s2);
    v2 = v(i);
    s2 = s(i);
end

[v, s] = TFitting_2RVs(v1,v2,s1,s2);

%%
%[ref] Fitting the Distribution of Linear Combinations of $t-$Variables
%with more than 2 Degrees of Freedom (submitted to Journal of Probabilities and Statistics)
end