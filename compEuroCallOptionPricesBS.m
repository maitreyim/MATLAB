%Problem 3, Part B
% Black-Scholes Pricing Model of European Call Option
% no of ietration
% c, % call price
% r, % rate of interest
% si, % standard deviation
% S0, % initial stock price
% T,  % time to expiry
% K, % strike price )
function [ callVal ] = compEuroCallOptionPricesBS(S0,T, X, si, r) 
d1 = 0.0498673470;
d2 = 0.0211410061; 
d3 = 0.0032776263;
d4 = 0.0000380036;
d5 = 0.0000488906;
d6 = 0.0000053830;

dd1 = (log(S0/X) + (r+(si^2)/2)*T)/si*sqrt(T);
dd2 = dd1 - si*sqrt(T);

if dd1 >= 0 
    Ndd1 = 1-(1/2)*(1+d1*dd1 + d2*dd1^2 + d3*dd1^3 + d4*dd1^4 + d5*dd1^5 + d6*dd1^6)^(-16);
else 
    dd1 = (-1) * dd1;
    Ndd1 = (1/2)*(1+d1*dd1 + d2*dd1^2 + d3*dd1^3 + d4*dd1^4 + d5*dd1^5 + d6*dd1^6)^(-16);
end

if dd2 >= 0
    
    Ndd2 = 1-(1/2)*(1+d1*dd2 + d2*dd2^2 + d3*dd2^3 + d4*dd2^4 + d5*dd2^5 + d6*dd2^6)^(-16);
else 
    dd2 = (-1) * dd2;
    Ndd2 = (1/2)*(1+d1*dd2 + d2*dd2^2 + d3*dd2^3 + d4*dd2^4 + d5*dd2^5 + d6*dd2^6)^(-16);
end

c = S0*Ndd1 - X*exp(-r*T)*Ndd2;
callVal = c;    
end

