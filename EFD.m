function [Prices, St, gridS] = EFD( type, callput, S0, K, r, sigma, T, delta, delta_t, Smax, Smin, logPrice )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function called EFD is written to evaluate option pricing using
%Explicit Finite Difference Method.The Inputs are as follows:
%Type American (A) or European (E)
%callput is Call (C) or Put (P)
%S0 is stock price, K is strike price, r is interest rate, sigma is
%volatility, T is time to maturity
%Delta is the increment in stock price based on S or X = ln(S)
%Delta_t is the time step size in the grid
%Smax and Smin are the maximum and minimum of the underlying prices we want
%to find the option prices for.
%LogPrices equal to 1 if we are going to use X = ln(S) to generate the stock
%prices else equal to 0 if using S.
%OUTPUT of the Function EFD:
%OptionPrices is the prices of all the options in the grid at each time
%step for different underlying price.
%St is the prices of the underlying security
%gridS is the row in the grid where the current stock price is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
columns = ceil(T/delta_t);%Time steps
%CASE I:If we are going to use X = ln(S) or S to generate stock prices
if(logPrice == 1)
    [Xt, rows, gridS] = StockDynamics_FD(log(S0), log(Smin), log(Smax), delta, 0);
  %CASE II: We are going to use ordinary stock price
    St = exp(Xt);%Change the prices from log
else
    [St, rows, gridS] = StockDynamics_FD(S0, Smin, Smax, delta, 0);
end;
opt = zeros(rows, columns); %Create array to capture option prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERTION OF THE GRID/MATRIX
%Option prices at maturity
if(callput == 'C')
   opt(:,columns) = max(St - K, 0);%%%%% FOR CALL
else
   opt(:,columns) = max(K - St, 0); %%%%FOR PUT
end;
%%%%%CASE I:
%%%EXPLICIT FINITE DIFFERENCE PARAMETERS WITH X=Ln(S)
%Explicit X, Pu, Pm and Pd based on when prices are generated with X = ln(S)
Pu = delta_t*(sigma^2/(2*delta^2)+(r-sigma^2/2)/(2*delta));
Pm = 1-delta_t*sigma^2/(delta^2); %-r*delta_t;
Pd = delta_t*(sigma^2/(2*delta^2)-(r-sigma^2/2)/(2*delta));
%Here we calculate the price in each time step and node backward.
%We start from the second last column where we have calculated the price in 
%last column above.
for i = columns-1:-1:1
    EV = 0; %Exercise value
    %Go through every row in each column
    for j = 1:rows
%%%%%%% CASE II:
%If we didn't use log prices then we need to have other Pu, Pm and
        %Pd
        if(logPrice == 0)
%%%EXPLICIT FINITE DIFFERENCE PARAMETERS WITH ORIGINAL B-S PDE
            Pu = delta_t*((r*j)/2+(sigma^2*j^2)/2);
            Pm = 1-delta_t*sigma^2*j^2;
            Pd = delta_t*((-r*j)/2+(sigma^2*j^2)/2);
        end

        if(j == 1) %If first row in each column, we need to set bondaries
            if(callput == 'C') 
                opt(j,i) = St(j)-K*exp(-r*(columns-i)*delta_t);
                EV = St(j) - K;
            else
               opt(j,i) = 0;
                EV = K - St(j);
            end;
        elseif(j == rows) %Other bondaries for the last row in each column
            if(callput == 'C') 
               opt(j,i) = 0;
                EV = St(j) - K;
            else
                opt(j,i) = K*exp(-r*(columns-i)*delta_t);
                EV = K - St(j);
            end;
        else %If not bondaries we calculate according to the explicit formula
            if(callput == 'C'); EV = St(j) - K;
            else EV = K - St(j); end;
           opt(j,i) = (1/(1+r*delta_t))*(Pu*opt(j-1,i+1)+Pm*opt(j,i+1)+Pd*opt(j+1,i+1));
        end;
        %THE AMERICAN CONDITION
        if(type == 'A' && EV > opt(j,i))
            opt(j,i) = EV;
        end;
    end;
end;
Prices = opt; %return option prices in each node.
end

