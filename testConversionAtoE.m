function [X] = EurpeanOptionOnCommodity(isCall, S, X, r, b, T, sigma)
%--------------------------------------------------------------------------
% Calculates price of european option on a commodity
% 
% isCall        = Call = 1, Put = 0
% S             = Price of underlying asset (S=F option on future)
% X             = Strike Price of Option
% r             = Risk free interest rate
% b             = Cost-of-carry rate (b=0 for option on future)
% T             = Time to Maturity
% sigma         = Volatility
%--------------------------------------------------------------------------

d1 = (log(S/X) + (b+0.5*sigma^2)*T) / (sigma * sqrt(T));
d2 = d1 - sigma *sqrt(T);

% Value of european call option on commodity
if (isCall)
    optionPrice = S * exp((b-r)*T) * normcdf(d1,0,1) - X * exp(-r*T) * normcdf(d2,0,1);
else
    optionPrice = X * exp(-r*T) * normcdf(-d2,0,1) - S * exp((b-r)*T) * normcdf(-d1,0,1);
end

X = optionPrice;