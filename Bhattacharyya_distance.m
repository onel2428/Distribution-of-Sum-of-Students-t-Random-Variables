function divBha = Bhattacharyya_distance(t,s,v)

%% Inputs
% t      -> samples from the sum of K t-RVs
% s      -> vector with the scales of the K t-RVs
% v      -> vector with the number of degrees of freedom of the K t-RVs

%% Outputs
% divBha -> Bhattacharyya distance metric


%% The computation of the Bhattacharyya distance metric is based on the procedure described in Section 4 in [ref]
[N,edges] = histcounts(t1,'normalization','pdf');

edges(1)=[];
y = @(x) pdf('tLocationScale',x,0,s,v);

for i=1:length(t)
    [~,p]=min(abs(t(i)-edges));
    f(i) = N(p);
end

f(f==0)=1e-50;
f = f/sum(f);
yt = y(t);
yt(yt==0)=1e-50;
yt = yt/sum(yt);

divBha = -log(sum(sqrt(f.*yt)));

%%
%[ref] Fitting the Distribution of Linear Combinations of $t-$Variables
%with more than 2 Degrees of Freedom (submitted to Journal of Probabilities and Statistics)

end