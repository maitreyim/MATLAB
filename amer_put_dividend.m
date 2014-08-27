%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%American Put with dividend%%%%%
function V = amer_put_dividend(S,E,r,sigma,q,tau,td,M)
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
payoff = max(E - Si,0);
Vi = max(payoff,Vi);
if i == m
Si = Si/(1-q);
payoff = max(E - Si,0);
Vi = max(payoff,Vi);
end
end
V = Vi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AMerican call with Fix Cash dividend
function [C0] = AmFixDiv (S0 ,X,r,tau ,sigma ,D,tauD ,M)
% compute constants
 f7 = 1; dt = tau / M; v = exp (-r * dt);
u = exp ( sigma * sqrt (dt)); d = 1 /u;
p = ( exp (r * dt) - d) / (u - d);
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
C(f7+j) = max (C(f7 + j), S(f7+j) - X);
else
C(f7+j) = max (C(f7 + j), S(f7+j) + D* exp (-r * (tauD -t)) - X);
end
end
end
C0 = C(f7 +0) ;
%%%%%%%%%%%%%%%%%%%%%%%%%
%Am Call with Div Yield
function [Binomial] = BinomialAmerican(CallPut, AssetP, Strike, RiskFree, Div, Time, Vol, nSteps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Cox, Ross & Rubinstein (1979) Binomial Tree for American Call/Put Option Values based
% on the following inputs:
% CallPut           =       Call = 1, Put = 0
% AssetP            =       Underlying Asset Price
% Strike            =       Strike Price of Option
% RiskFree          =       Risk Free rate of interest
% Div               =       Dividend Yield of Underlying
% Time              =       Time to Maturity
% Vol               =       Volatility of the Underlying
% nSteps            =       Number of Time Steps for Binomial Tree to take
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = Time / nSteps;

if CallPut
    b = 1;
end
if ~CallPut
    b = -1;
end

RR = exp(RiskFree * dt);
Up = exp(Vol * sqrt(dt));
Down = 1 / Up;
P_up = (exp((RiskFree - Div) * dt) - Down) / (Up - Down);
P_down = 1 - P_up;
Df = exp(-RiskFree * dt);

% Sets up the asset movements on the binomial tree
for i = 0:nSteps
    State = i + 1;
    St = AssetP * Up ^ i * Down ^ (nSteps - i);
    Value(State) = max(0, b * (St - Strike));
end

% Works backwards recursively to determine the price of the option    
for TT = (nSteps - 1):-1:0
    for i = 0:TT
        State = i + 1;
        Value(State) = max((b * (AssetP * Up ^ i * Down ^ (abs(i - TT)) - Strike)), (P_up * Value(State + 1) + P_down * Value(State)) * Df);
    end
end

Binomial = Value(1)

 