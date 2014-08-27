% Question 1, Partd
clear all
S0    = 32;
T     = 0.5;
r     = .05;
sigma = .24;
K     = 30;
M = 10;
delta = T / M;
p = .5 + r * sqrt(delta)/(2*sigma);
assert(p < 1)
u=exp(sigma*sqrt(delta));
d=exp(sigma*sqrt(delta))*(-1);
p= 0.5+0.5*(((r-(sigma^2)/2)*sqrt(delta))/sigma);
V = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
for m = 1:M
    V = 1/(1 + r * delta ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );
end
V