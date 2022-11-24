function [vf, sf, rf] = TFitting_KRVs_CharacteristicFunction(v,s,r)

%% Inputs
% v   -> vector with the degrees of freedom of each of the input t-RVs
% s   -> vector with the scale of each of the input t-RVs
% r   -> parameter of the characteristic function to be used for the fitting.
% (Different parameters r may lead to different fitting accuracy figures)

%% Outputs
% vf  -> degree of freedom of the fitting distribution 
% sf  -> scale of the fitting distribution
% rf  -> parameter r of the characteristic function according to eq. (40) in
% [ref]


%% Relevant Statistics for the Fitting
K = length(v);
M2 = SecondMoment(v,s);  % Second Moment
CF_func = @(rx) 2^K*(rx/2).^(sum(v)/2).*prod((sqrt(v).*s).^(v/2).*besselk(v/2,sqrt(v).*s.*rx)./gamma(v/2)); % Characteristic function (CF)
rf = M2^(-1/2);          % Parameter r proposed in eq. (40) in [ref]

CF(1) = CF_func(r);  % CF value for the input r parameter
CF(2) = CF_func(rf); % CF value for the r proposed in eq. (40) in [ref]

%% Fitting based on the CF and second moment for both the parameter r provided as input and that proposed in eq. (40) in [ref]
% Check the procedure in Section 3.2 in [ref]
fv = @(vx)  (r*sqrt((vx-2)*M2)).^(vx/2).*besselk(vx/2,r*sqrt((vx-2)*M2))./(2.^(vx/2-1).*gamma(vx/2))-CF(1);

vf(1) = fzero(fv,[2.00001,160]);
sf(1) = sqrt((vf(1)-2)*M2/vf(1));

fvf = @(vx)  (rf*sqrt((vx-2)*M2)).^(vx/2).*besselk(vx/2,rf*sqrt((vx-2)*M2))./(2.^(vx/2-1).*gamma(vx/2))-CF(2);
vf(2) = fzero(fvf,[2.00001,160]);
sf(2) = sqrt((vf(2)-2)*M2/vf(2));

%%
%[ref] Fitting the Distribution of Linear Combinations of $t-$Variables
%with more than 2 Degrees of Freedom (submitted to Journal of Probabilities and Statistics)
end