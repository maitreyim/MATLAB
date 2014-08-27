% Question 4
% Estimation of European and American Put Pricing
% European Put Pricing
S0    = 80;
T     = 1;
r     = .05;
sigma = .03; 
K=100;
M = 20;
delta =T/M;
p = .5 + r * sqrt(delta)/(2*sigma);
u = 1 + sigma * sqrt(delta);
d = 1 - sigma * sqrt(delta);

C = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
for m = 1:M
    C = 1/(1 + r * delta ) * ( C(1:(end-1)) * p + C(2:end) * (1-p) );
end
C
%plot
% American Put Pricing
RR = exp(r  * delta);
Up = exp(sigma* sqrt(delta));
Down = 1 / Up;
P_up = (exp((r ) * delta) - Down) / (Up - Down);
P_down = 1 - P_up;
Df = exp(-r  * delta);
% Sets up the asset movements on the binomial tree
for i = 0:M
    State = i + 1;
    St = S0 * Up ^ i * Down ^ (M - i);
    Value(State) = max(0, (-1) * (St - K));
end
% Works backwards recursively to determine the price of the option    
for TT = (M - 1):-1:0
    for i = 0:TT
        State = i + 1;
        Value(State) = max(((-1) * (S0 * Up ^ i * Down ^ (abs(i - TT)) - K)), (P_up * Value(State + 1) + P_down * Value(State)) * Df);
    end
end
%plot


