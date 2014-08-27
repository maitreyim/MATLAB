%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Computational Finance Project # 6
%%%Finite Difference Method
%%%%%Submitted by Maitreyi Mandal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Initialization as per the Given Problem
r = 0.04; %Interest rate
sigma = 0.2; %Volatility
S0 = 10; %Current Stock Price, the middle price 
K = 10; %Strike Price
delta_t = 0.002; %Delta t
T = 0.5; %Time to maturity
%Delta X is used for both S and and X. We can only us S for Explicit FDM
delta_x =sigma * sqrt(3*delta_t);
%If equal to 1 then we are going use ln(X) to generate the stock prices
%if equal to 0 then we use S. Only us S for Explicit FDM
IndicatorX = 0; 
%StockPrice Range and the corersponding increment
Smax = 20;
Smin = 1;
increments = 1;
type = 'E'; %A for American, E for European
callput = 'C'; %C for Call and P for Put
%The three Finite Difference Methods that we are using here, need to
%comment out the ones which we are not using:
%[option, St, gridS] = IFD(type, callput, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin );
[option, St, gridS] = EFD(type, callput, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin, IndicatorX);
%[option, St, gridS] = CrankNicolsonFiniteDifferenceMethod(type, callput, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin );
%To get how many prices we have
stockpriceArray = (Smax - Smin + 1) / increments;
%Table to compare prices from different methods
values = zeros(stockpriceArray, 2);
for i = stockpriceArray:-increments:1
    %Aim of this step is to minimize the error to get the price closest to the prices we are
    %looking for as an underlying stock price which should be between $1 to
    %$20 as given in the problem
    delta = 1000;
    for j = 1:size(St,1)%Prices that we want
        StPrice = St(j);%Prices that we have now
        if(abs(StPrice - i) < delta)
            delta = abs(StPrice - i);
            values (i, 1) = StPrice;%First column is the price of underlying
            values (i, 2) = option(j);%Second column is the option price from FDM
                    end;
    end;
end;
vals = values%Get the array of stock price and corresponding option prices as wanted in the problem.