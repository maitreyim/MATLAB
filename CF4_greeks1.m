function [ V1] = CF4_greeks1(sigma,delta,r )
delta = T / M;
S0=49;
epsillion=0.05;
S1=S0+epsillon;
u=exp(r*delta)*(1+sqrt(exp((sigma^2)*delta)-1));
d=exp(r*delta)*(1-sqrt(exp((sigma^2)*delta)-1));
p=0.5;

V1 = max(S1 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
for m = 1:M
    V = 1/(1 + r * delta ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );
end

V1