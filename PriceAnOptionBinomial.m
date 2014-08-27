function [ optionValue ] = PriceAnOptionBinomial(type, kind, S0, K, r, sigma, T, n, u, d, p)
%Funtion to value an option using binomial tree.
%Type American (A) or European (E)
%Kind is Call (C) or Put (P)
%S0 is stock price, K is strike price, r is interest rate, sigma is
%volatility, T is time to maturity, n number of steps in the tree, u and d
%are how much the stock will go up or down in each step, p is the likelihood that stock will go up
%   Detailed explanation goes here
x = n + 1; %Because we want to include S0 and V0 in the arrays we need one extra column
%Create array to accumulate stock prices and option prices
St = zeros(x,x);
Vt = zeros(x,x);
St(1,1) = S0; %First stock price is the current price
delta = T/n; %Step size
%Check if all the constrains are ok
if(S0 >= 0 && K >= 0 && sigma >= 0 && T >= 0 && n >= 0)
    %Fill the stock price array of all possible stock prices
    for column = 2:x
        uPower = column-1;
        dPower = 0;

        for row = 1:column
            St(row, column) = u^uPower * d^dPower * S0;
            uPower = uPower - 1; %How much price goes up
            dPower = dPower + 1; %How much price goes down
        end;
    end;
    %Then go backward and calculate the option price from the stock prices
    for column = x:-1:1
        for row = 1:column
            %If last column take the max of strike price and stock price
            if(column == x)
                if(kind == 'C') %If call
                    Vt(row,column) = max(St(row,column) - K, 0);
                elseif(kind == 'P') %If put
                    Vt(row,column) = max(K - St(row,column), 0);
                end;
            else %Else get the value in each node
                node1 = Vt(row,column+1);
                node2 = Vt(row+1,column+1);
                %And calculate the continuous value (CV)
                CV = exp(-r*delta)*(p * node1 + (1-p) * node2);
                %Calculate the exercise value (EV)
                if(kind == 'C')
                    EV = max(St(row,column) - K, 0);
                elseif(kind == 'P')
                    EV = max(K - St(row,column), 0);
                end;
                value = CV;
                %Only use EV if American option and EV > CV
                if(type == 'A' && CV < EV)
                    value = EV;
                end;

                Vt(row,column) = value;
            end;
        end;
    end;

    optionValue = Vt(1,1);
else
    optionValue = 0; %If constrain above not ok, value = 0
end;
end

