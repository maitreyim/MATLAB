function Main( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TERM PAPER FOR OPTIONS MARKET CLASS
%REPLICATING PAPER BY KEMNA & VORST ON AVERAGE VALUE OPTION
%CODE SUBMITTED BY: MAITREYI MANDAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PARAMETERS AS IN THE ORGINIAL PAPER
S0 = 40;%Stock Price
K = [35, 40, 45];%Exercise Rate
r = [0.03, 0.05, 0.07];%Risk-free interest rates
sigma = [0.2, 0.3, 0.4];%Volatilities as used in the paper
T0 = 0; %Initial Time
T = 4/12;
n = 88;
N = 10000;% 
PricesArrayOpt = nan(27,7);
count = 1;
for i=1:3
    for j=1:3
        for k=1:3
            [Ca, Cb, SEb, Cc, SEc] = ValueAverageOption( S0, K(k), r(i), sigma(j), T0, T, n, N );
            priceArrayOpt(count,:) = [r(i),sigma(j),K(k),Ca,Cb,SEb,Cc,SEc];
            count = count + 1;
        end
    end
end
end

function [ Ca, Cb, SEb, Cc, SEc ] = ValueAverageOption( S0, K, r, sigma, T0, T, n, N )
%ASIAN/AVERAGE OPTION VALUATION USING ARITHMATIC & GEOMETRIC MEAN
%   
Ave = nan(N,1);%Pre-allocating memory space by defining a matrix of zeros
Geo = nan(N,1);%Pre-allocating memory space by defining a matrix of zeros
Z = randn(n,N);% Generating random number to be used in the monte carlo simulation

for i = 1:N
    St = StockPrices( S0, r, sigma, T0, T, n, Z(:,i) );
    Ave(i) = exp(-r*(T-T0))*max(mean(St)-K,0);
    Geo(i) = AverageOptionVarianceReduction(S0, St, K, r, sigma, T0, T, Ave(i));
end;

Ca = blsprice(S0, K, r, T, sigma);%Ca gives ordinary European Call Prices using Black-Scholes
Cb = mean(Ave);%Cb gives Asian/Average Prices using Arithmetic Average
SEb = std(Ave)/sqrt(N);%SEb gives Standard Errors/Deviation in Cb
Cc = mean(Geo);%
SEc = std(Geo)/sqrt(N);
end

function [ St ] = StockPrices( S0, r, sigma, T0, T, n, Z )
%GENERATION OF STOCK PRICES, PAPER USES A LOG PROCESS BUT HERE A MORE
%GENERAL APPROACH WAS TAKEN
St = zeros(n+1,1);
St(1) = S0;
for i = 2:n+1
    St(i,1) = St(i-1,1)*exp((r-0.5*sigma^2)*(T-T0)/n+sigma*sqrt((T-T0)/n)*Z(i-1));
end;

end

function [ valueVarRed ] = AverageOptionVarianceReduction( S0, St, K, r, sigma, T0, T, Y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d1 = (log(S0/K)+0.5*(r+(sigma^2)/6)*(T-T0))/(sigma*sqrt((T-T0)/3));
d2 = d1-sigma*sqrt((T-T0)/3);
dStar = 0.5*(r-(sigma^2)/6)*(T-T0);

EW = exp(dStar)*S0*normcdf(d1,0,1)-K*normcdf(d2,0,1);
W = exp(-r*(T-T0))*max(geomean(St)-K,0);

valueVarRed = exp(-r*(T-T0))*EW + (Y-W);
end
