function M = SecondMoment(nu,sigma)

%% Inputs
% nu    -> vector with the number of degrees of freedom of the t-RVs
% sigma -> vector with the scale of the t-RVs

%% Output
% M     -> second moment of the sum of the t-RVs

% Second moment computation
M = sum(sigma.^2.*nu./(nu-2));


end