function V = amer_call_dividend(S,E,r,sigma,q,tau,td,M)
dt = tau/M;
A = (exp(-r*dt) + exp((r+sigma^2)*dt))/2;
u = A + sqrt(A^2-1);
d = 1/u;
p = (exp(r*dt)-d) / (u-d);
Si = S*(1-q)*(u.^(M:-2:-M));
Vi = max(E - Si,0);
m = find((tau-dt:-dt:0) == td);
for i = 1:M
Vi = exp(-r*dt)*(p*Vi(1:end-1)+(1-p)*Vi(2:end));
Si = Si(1:end-1)/u;
payoff = max(0,Si-E);
Vi = max(payoff,Vi);
if i == m
Si = Si/(1-q);
payoff = max(E - Si,0);
Vi = max(payoff,Vi);
end
end
V = Vi;
end

 