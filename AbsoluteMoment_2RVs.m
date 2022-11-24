function M = AbsoluteMoment_2RVs(nu,sigma)

%% Inputs
% nu    -> two-element vector with the number of degrees of freedom of two t-RVs
% sigma -> two-element vector with the scale of two t-RVs

%% Output
% M     -> absolute moment of the sum of the two t-RVs

%% The absolute moment computation is based on the procedure described in Section 3.1 in [ref]
v1=nu(1);
v2=nu(2);

s1=sigma(1);
s2=sigma(2);

alf = @(v) gamma((v+1)/2)/(gamma(v/2)*sqrt(v*pi));

I1 = alf(v1)*alf(v2)*beta(1/2,(v1+v2-1)/2)*v1^(3/2)*s1*...
    hypergeom([(v2+1)/2 1/2],(v1+v2)/2,1-s1^2*v1/(s2^2*v2)) ...
    /((v1-1)*s2);
Ia1 = sqrt(v2)*gamma((v2-1)/2)/(2*sqrt(pi)*gamma(v2/2));
w1 = v2^((v2+1)/2)*(v1*s1^2/s2^2)^((1-v2)/2);
w2 = v2*s2^2/(v1*s1^2)-1;
Ia2=Ia1-0.5*alf(v2)*w1*integral(@(z) z.^(-3).*betainc(z.^2,v1/2,1/2).*(w2+1./z.^2).^(-(v2+1)/2),0,1);

I22p = w1*integral(@(z) z.^(-3).*betainc(z.^2,v1/2,1/2).*(w2+1./z.^2).^(-(v2+1)/2),0,1);
I22pv1v2 = v1*(1-sqrt(pi)*gamma(v1-1/2)*2^(1-v1)/gamma(v1/2)^2)/(v1-1);


I2 = 2*Ia2-Ia1;

M=2*(s1*I1+s2*I2);
%Me = integral2(@(x1,x2) abs(s1*x1+s2*x2).*pdf('tLocationScale',x1,0,1,v1).*pdf('tLocationScale',x2,0,1,v2),-Inf,Inf,-Inf,Inf);

%%
%[ref] Fitting the Distribution of Linear Combinations of $t-$Variables
%with more than 2 Degrees of Freedom (submitted to Journal of Probabilities and Statistics)

end