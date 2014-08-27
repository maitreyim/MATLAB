%Question 3
K=50;
r=0.03;
sigma=0.2;
T=0.3846;
M = 10;
delta = T / M;
p = .5 + r * sqrt(delta)/(2*sigma);
assert(p < 1)
c=0.5*(exp(-r*delta))+exp(r+(sigma^2)*delta);
d=c-sqrt(c^2-1);
u=1/d;
p=((exp(r*delta)-d)/(u-d));
S0=10;
increment = 0.0005;
delta = zeros(1, 40);
gamma = zeros(1, 40);
theta = zeros(1, 40);
vega = zeros(1, 40);
rho = zeros(1, 40);
for m = 1:M
for i = 15:1:25
    S0 = i;
    %Computing Delta: the First Option Greek, which measures sensitivity of
    %change in call price w.r.t change in stock price, this is equal to
    %N(d1)
    S1 = S0;
    S2 = S0 + increment;
    V = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
    
    V = 1/(1 + r * delta ) * ( V(1:(end-1)) * p + V(2:end) * (1-p) );
    
    V;
%     First_S1 = blsprice(S1, X, r, T, sigma);
%     Next_S2 = blsprice(S2, X, r, T, sigma); 
    delta1 = (Next_S2 - First_S1) / increment;
    delta(i - 14) = delta1;
end;