function [C0] = AmericanCallFixDiv (CallPut,S0 ,X,r,tau ,sigma ,D,tauD ,M)
% compute constants
 f7 = 1; dt = tau / M; v = exp (-r * dt);
u = exp ( sigma * sqrt (dt)); d = 1 /u;
p = ( exp (r * dt) - d) / (u - d);
if CallPut
    b = 1;
end
if ~CallPut
    b = -1;
end
% adjust spot for dividend
S0 = S0 - D * exp (-r * tauD );
% initialise asset prices at maturity ( period M)
S = zeros (M + 1 ,1);
S(f7 +0) = S0 * d^M;
 for j = 1:M
S(f7+j) = S(f7+j - 1) * u / d;
end
% initialize option values at maturity ( period M)
C = max (S - X, 0);
% step back through the tree
for i = M -1: -1:0
for j = 0:i
C(f7+j) = v * ( p * C(f7+j + 1) + (1-p) * C(f7+j));
S(f7+j) = S(f7+j) / d;
t = tau * i / M;
if t > tauD
C(f7+j) = max (C(f7 + j), S(f7+j) - b*X);
else
C(f7+j) = max (C(f7 + j), S(f7+j) + D* exp (-r * (tauD -t)) - b*X);
end
end
end
C0 = C(f7 +0) ;
