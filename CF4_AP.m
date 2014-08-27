
function [AP] = CF4_AP(S0,T,r,sigma,K )
    
% Question 4
% American Put Pricing
M = 150;
delta =T/M;
p = .5 + r * sqrt(delta)/(2*sigma);
u = 1 + sigma * sqrt(delta);
d = 1 - sigma * sqrt(delta);
AP = max(S0 * (-1)*( (d.^(0:M)) .* (u.^(M:(-1):0)) )' +K, 0);
for m = 1:M
    CV = 1/(1 + r * delta ) * ( AP(1:(end-1)) * p + AP(2:end) * (1-p) );
    EV = max( K - S0.*(d.^(0:(M-m))).*(u.^((M-m):-1:0)),0 )';
    AP = max(CV, EV);
end
AP;
end