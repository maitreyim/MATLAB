S0    = 615;
T     = 0.5;
r     = .05;
sigma = .20;
N=7000;
M = 20;
delta = T / M;
% p = .5 + r * sqrt(delta)/(2*sigma);
% assert(p < 1)
c=0.5*(exp(-r*delta))+exp(r+(sigma^2)*delta);
d=c-sqrt(c^2-1);
u=1/d;
p=((exp(r*delta)-d)/(u-d));
V=ones(1,7000);
K=ones(1,7000);
for i=20:N
    K(i)=i;
% Calculate Binomial Call at each node
V(i) = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K(i), 0);
%Backward Recursion
for m = 1:M
    V = 1/(1 + r * dt ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );
end

V
