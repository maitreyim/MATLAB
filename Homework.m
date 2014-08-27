%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Homework - Derivatives
%05/16/2012
%Code Written by:Maitreyi Mandal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partameters as used in Geske/Shastri 1985 paper

S0 = 40;% Stock Price
sigma = 0.3;% volatility
r = 0.05;% Risk Free Rate
K = [35, 40, 45]; % Exercise or Strike Price
Time = [1/12, 4/12, 7/12];% Time to Expiry
%We need to choose between dividend or dividendYield. Either one needs
%at least to be equal to one. Both zero if no dividents. 
divident = 0.50;% Fixed Cash Dividend or BFCD
dividentYield = 0; % Fixed Dividend Yield or BFDY
%dTime needs to be zero if we have no divident not dividentYield
dTime = [0.5/12, 3.5/12, 6.5/12];% Time when Dividends are paid
n =140;
optionValue = zeros(3,3);
for i = 1:size(Time, 2)
    T = Time(i);
    delta = T / n;
    d = exp(-sigma*sqrt(delta));
    u = exp(sigma*sqrt(delta));
    p = 1/2+ 1/2*(((r-sigma^2/2)*sqrt(delta))/sigma);
    for j = 1:size(K,2)
        Kt = K(j);
        optionValue(j,i) = PriceAnOptionBinomial('A', 'C', S0, Kt, r, sigma, T, n, u, d, p, dTime, divident, dividentYield);        
    end;
end;

put = optionValue
