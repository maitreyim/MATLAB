function [optionPrice, optionVola] = convertAmericanToEuropean(CallOrPut, optionPrice, S, X, r, T)
%--------------------------------------------------------------------------
% Converts a set of american options into european options using 
% Barone-Adesi and Whaley (1987) to get implied vola and from there 
% calculate corresponding euroean options with Black formula
% -> See "Unspanned Stochastic Volatility [..]", Trolle, Schwartz (2009)
%
% Function returns european option and implied volatility
%
% CallOrPut     = Call = 1, Put = 0
% optionPrice   = Market price of american option
% S             = Price of underlying asset (S=F option on future)
% X             = Strike Price of Option
% r             = Risk free interest rate
% T             = Time to Maturity
%
% Return Values
% optionPrice   = European call/put option price
% optionVola    = Vola of option
%--------------------------------------------------------------------------


    try
        impliedVol = impliedVolatility_BAW(CallOrPut, optionPrice, S, X, r, 0.0, T);
        europeanOption = EurpeanOptionOnCommodity(CallOrPut, S, X, r, 0.0, T, impliedVol);
    catch
        impliedVol = NaN;
        europeanOption = NaN;
    end

    optionPrice = europeanOption;
    optionVola = impliedVol;

end

