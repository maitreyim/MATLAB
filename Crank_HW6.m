function [ optionPrices, St, middle] = Crank_HW6( type, kind, S0, K, r, sigma, T, delta_x, delta_t, Smax, Smin )
%Funtion to valuate an option using Crank-Nicolson Finite Difference Method.
%INPUT
%Type American (A) or European (E)
%Kind is Call (C) or Put (P)
%S0 is stock price, K is strike price, r is interest rate, sigma is
%volatility, T is time to maturity
%Delta_x is the increment in stock price based on X = ln(S)
%Delta_t is the time step size in the grid
%Smax and Smin are the maximum and minimum of the underlying prices we want
%to find the option prices for.
%OUTPUT
%OptionPrices is the prices of all the options in the grid at each time
%step for different underlying price.
%St is the prices of the underlying security
%Middle is the row in the grid where the current stock price is

columns = ceil(T/delta_t);%How many columns/steps we need.
%Get log stock prices
[Xt, rows, middle] = StockPricesFiniteDifferenceMethod(log(S0), log(Smin), log(Smax), delta_x, 0);
St = exp(Xt);%Change the prices from log
%Matrices we need to use to get the option prices for every row/node at
%each time step.
A = zeros(rows, rows);
o = zeros(rows, 1);%Option prices at time t
b = zeros(rows, 1);%Option prices at time t + 1

%Crank-Nicolson X, Pu, Pm and Pd based on when prices are generated with X = ln(S)
v = r - sigma^2/2;
Pu = (-0.25)*delta_t*(sigma^2/delta_x^2+v/delta_x);
Pm = 1+delta_t*sigma^2/(2*delta_x^2)+r*delta_t/2;
Pd = (-0.25)*delta_t*(sigma^2/delta_x^2-v/delta_x);

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
if(kind == 'C')
    o = max(St - K, 0);
else
    o = max(K - St, 0);
end;
%Here we calculate the price in each time step and node backward.
for i = columns:-1:1
    EV = 0;%Exercise value
    b = o;%Option value at time t+1
    if(i ~= columns)%If last column(maturity) then skip. Caluclated above
        for j = 1:rows
            %Set bondaries based on Call(C) or Put(P)
            if(kind == 'C')
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
                b(j) = -Pu*o(j-1)-(Pm-2)*o(j)-Pd*o(j+1);
            end;
            %If we have an American option
            if(type == 'A' && EV > b(j)); b(j) = EV; end;
        end;
    end;
    
    o = A\b;%Get option prices at time t
end;

optionPrices = o;%Return option prices

end



