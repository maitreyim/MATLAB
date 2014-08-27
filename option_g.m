function [V] = option_g(S0,sigma,T,r,K )
%Function to be referenced in main script to generate the Option Greeks
M = 500;% Time Steps
delta = T / M;
u=exp(r*delta)*(1+sqrt(exp((sigma^2)*delta)-1));% Up movement
d=exp(r*delta)*(1-sqrt(exp((sigma^2)*delta)-1));%Down Movement
p=0.5;% probability
V = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);% Call Price at each node
for m = 1:M
    V = 1/(1 + r * delta ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );%Call Price at the First Node
end

V;





