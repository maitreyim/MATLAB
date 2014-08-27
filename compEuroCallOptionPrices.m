
% %no of ietration
%                                                       c, % call price
%                                                       r, % rate of interest
%                                                       si, % standard deviation
%                                                       S0, % initial stock price
%                                                       T,  % time to expiry
%                                                       K) { % strike price )

function [ callVal, callValAnti ] = compEuroCallOptionPrices(  n,r, si, S0, T,K)  
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

W = sqrt(T)*randn(n,1);
%Y = -X;
ST = S0*exp(si*W+ (r - (si^2.0)/2.0)*T);
c = exp(-r*T)*max(ST-K,0) ;
callVal = mean(c);    

STA = S0*exp(si*(-1)*W+ (r - (si^2.0)/2.0)*T);
cA = exp(-r*T)*max(STA-K,0) ;
cAvg = (c+cA)/2;
callValAnti = mean(cAvg);
end

