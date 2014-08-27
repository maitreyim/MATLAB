function [ Prices, St, gridS] = CrankNicolson( type, callput, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The Function CrankNicolson is written to evaluate the option price using
%Crank-Nicolson Method of Finite Difference, the various INPUTS and OUTPUTS
%are given below:
%INPUT of this finction are:
%Type American (A) or European (E)
%callput is Call (C) or Put (P)
%S0 is stock price, K is strike price, r is interest rate, sigma is
%volatility, T is time to maturity
%Delta_x is the increment in stock price based on X = ln(S)
%Delta_t is the time step size in the grid
%Smax and Smin are the maximum and minimum of the underlying prices we want
%to find the option prices for.
%OUTPUT of this function is:
%OptionPrices is the prices of all the options in the grid at each time
%step for different underlying price.
%St is the prices of the underlying security
%gridS is the row in the grid where the current stock price is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
columns = ceil(T/delta_t);%How many columns/steps we need.
%Get log stock prices
[Xt, rows, gridS] = StockDynamics_FD(log(S0), log(Smin), log(Smax), delta_x, 0);
St = exp(Xt);%Change the prices from log
%Matrices we need to use to get the option prices for every row/node at
%each time step.
A = zeros(rows, rows);
opt = zeros(rows, 1);%Option prices at time t
b = zeros(rows, 1);%Option prices at time t + 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Crank-Nicolson PARAMETERS
%Crank-Nicolson X, Pu, Pm and Pd based on when prices are generated with X = ln(S)
v = r - sigma^2/2;% Nu as found in the transformed Black-Scholes PDE
Pu = (-0.25)*delta_t*(sigma^2/delta_x^2+v/delta_x);
Pm = 1+delta_t*sigma^2/(2*delta_x^2)+r*delta_t/2;
Pd = (-0.25)*delta_t*(sigma^2/delta_x^2-v/delta_x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENERATION OF THE GRID and appropriate Boundary
%Conditions
%Fill the A array according to the notes for IFD
for j = 1:rows
    if(j == 1)
        A(j,j) = 1; A(j,j+1) = -1;
    elseif (j == rows)
        A(j,j-1) = 1; A(j,j) = -1;
    else
        A(j,j-1) = Pu; A(j,j) = Pm; A(j,j+1) = Pd;
    end;
end;
%Option prices at maturity
if(callput == 'C')
    opt= max(St - K, 0);
else
    opt= max(K - St, 0);
end;
%Here we calculate the price in each time step and node backward.
for i = columns:-1:1
    EV = 0;%Exercise value
    b = opt;%Option value at time t+1
    if(i ~= columns)%If last column(maturity) then skip. Caluclated above
        for j = 1:rows
            %Set bondaries based on Call(C) or Put(P)
            if(callput == 'C')
                if(j == 1)
                    b(1) = St(1) - St(2);
                elseif(j == rows)
                    b(rows) = 0;
                end;
                EV = St(j) - K;
            else
                if(j == 1)
                   b(1) = 0;
                elseif(j == rows)
                    b(rows) = -(St(rows) - St(rows-1));
                end;
                EV = K - St(j);           
            end;
            %If not bondaries we calculate according to the Crank-Nicolson formula
            if(j ~= 1 && j ~= rows)
                b(j) = -Pu*opt(j-1)-(Pm-2)*opt(j)-Pd*opt(j+1);
            end;
            %AMERICAN CONDITION
            if(type == 'A' && EV > b(j)); b(j) = EV; end;
        end;
    end;
    
    opt= A\b;%Get option prices at time t
end;
Prices = opt;%Return option prices

end



