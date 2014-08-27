% %Question 1, Part C
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
u=exp(r-(sigma^2)*delta)+sigma*sqrt(delta);
d=exp(r-(sigma^2)*delta)-sigma*sqrt(delta);
p=0.5;
V = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
for m = 1:M
    V = 1/(1 + r * delta ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );
end
V
