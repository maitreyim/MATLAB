function [ optionValue ] = PriceAnOptionTrinomial( type, kind, S0, K, r, sigma, T, n, u, d, pu, pd, logPrice)
%Funtion to value an option using trinomial tree.
%Type American (A) or European (E)
%Kind is Call (C) or Put (P)
%S0 is stock price, K is strike price, r is interest rate, sigma is
%volatility, T is time to maturity, n number of steps in the tree, u and d
%are how much the stock will go up or down in each step, pu and pd are the
%likelihood that stock will go up or down
%LogPrice should be equal to 'L' if we are using log prices in
%calculations
%   Detailed explanation goes here
pm = 1-pu-pd; %Likelihood that the stock will not move
rowsTotal = 2*n + 1; %Total of rows in the array
columnsTotal = n + 1; %Total of columns in the array
%Create array to accumulate stock prices and option prices
St = zeros(rowsTotal,columnsTotal);
Vt = zeros(rowsTotal,columnsTotal);
St(1,1) = S0; %First stock price is the current price
delta = T/n; %Step size
%Check if all the constrains are ok
if(S0 >= 0 && K >= 0 && sigma >= 0 && T >= 0 && n >= 0)
    %Fill the stock price array of all possible stock prices
    for column = 2:columnsTotal
        uPower = column-1;
        dPower = 0;
        rows = 2*(column-1) + 1;
        for row = 1:rows
            St(row, column) = u^uPower * d^dPower * S0;
            if(logPrice == 'L')%If using log prices
                St(row, column) = S0 + u*uPower + d*dPower;
            end;
                
            if(uPower ~= 0)
                uPower = uPower - 1; %How much price goes up
            else %OR
                dPower = dPower + 1; %How much price goes down
            end;
        end;
    end;
    %Then go backward and calculate the option price from the stock prices
    for column = columnsTotal:-1:1
        rows = 2*(column - 1) + 1;
        for row = 1:rows
            %If last column take the max of strike price and stock price
            if(column == columnsTotal) 
                if(kind == 'C') %If call
                    Vt(row,column) = max(St(row,column) - K, 0);
                elseif(kind == 'P') %If put
                    Vt(row,column) = max(K - St(row,column), 0);
                end;
            else %Else get the value in each node
                node1 = Vt(row,column+1);
                node2 = Vt(row+1,column+1);
                node3 = Vt(row+2,column+1);
                %And calculate the continuous value (CV)
                CV = exp(-r*delta)*(pu * node1 + pm * node2 + pd * node3);
                if(logPrice == 'L') %If using log prices
                    CV = -r*delta + pu * node1 + pm * node2 + pd * node3;
                end;
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


