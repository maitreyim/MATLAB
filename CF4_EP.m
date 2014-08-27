function [EP] = CF4_EP(S0,T,r,sigma,K )
% Question 4
% Estimation of European and American Put Pricing
% European Put Pricing
M = 500;% No of Time Steps
delta =T/M;
p = .5 + r * sqrt(delta)/(2*sigma);%
u = 1 + sigma * sqrt(delta);
d = 1 - sigma * sqrt(delta);
EP = max(S0 * (-1)*( (d.^(0:M)) .* (u.^(M:(-1):0)) )' +K, 0);
for m = 1:M
    EP = 1/(1 + r * delta ) * ( EP(1:(end-1)) * p + EP(2:end) * (1-p) );
end
EP;
end



