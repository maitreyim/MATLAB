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
delta_x =sigma * sqrt(4*delta_t);
%If equal to 1 then we are going use ln(X) to generate the stock prices
%if equal to 0 then we use S. Only us S for Explicit FDM
useX = 1; 
%StockPrice Range and the corersponding increment
Smax = 20;
Smin = 1;
increments = 1;
type = 'A'; %A for American, E for European
kind = 'P'; %C for Call and P for Put

%We need to comment out the method we don't want to use
%[option, St, middle] = IFD_HW6(type, kind, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin );
[option, St, middle] = EFD_HW6(type, kind, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin, useX);
%[option, St, middle] = CrankNicolsonFiniteDifferenceMethod(type, kind, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin );

%To get how many prices we have
sizeOfArray = (Smax - Smin + 1) / increments;
%Table to compare prices from different methods
val = zeros(sizeOfArray, 2);

for i = sizeOfArray:-increments:1
    %Minimize the error to get the price closest to the prices we are
    %looking for as an underlying stock price
    error = 1000;
    for j = 1:size(St,1)%Prices that we want
        StPrice = St(j);%Prices that we have
        if(abs(StPrice - i) < error)
            error = abs(StPrice - i);
            val(i, 1) = StPrice;%First column is the price of underlying
            val(i, 2) = option(j);%Second column is the option price from FDM
                    end;
    end;
end;

option = val%Get the comparison matrix.