function [ Xt, rows, middle] = StockPricesFiniteDifferenceMethod(X0, Xmin, Xmax, delta_x, u)
%Function to calculate stock prices for use in Finite Difference Methods
%INPUT
%X0 is the central price
%Xmax and Xmin are the maximum and minimum of the underlying prices we want
%Delta_x is the increment in stock price
%Usin u is another way to get increment in stock prices
%OUTPUT
%Xt includes all the stock prices in the grid at each time
%rows is how many rows are in the grid
%Middle is the row in the grid where the current stock price is

%Get the range we need so we cover all prices(Xmax and Xmin) that we want.
range = max(X0-Xmin, Xmax-X0)/delta_x;
N = ceil(range); %Be sure to have integer
rows = 2*N + 1;%We will have N numbers on both side of the middle
middle = N + 1;

Xt = zeros(rows, 1);
Xt(middle, 1) = X0;%Set the current price

%Create stock prices
for i = 1:N
    if(u == 0)
        Xt(middle+i,1) = X0 - i*delta_x;
        Xt(middle-i,1) = X0 + i*delta_x;
    else
        Xt(middle+i,1) = X0*u^(-i);
        Xt(middle-i,1) = X0*u^(i);
    end;
end;

end

