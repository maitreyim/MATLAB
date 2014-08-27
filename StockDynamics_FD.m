function [ Xt, rows, gridS] = StockDynamics_FD(X0, Xmin, Xmax, delta_x, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function StockDynamics_FD is written to calculate stock prices for use in
%Finite Difference Methods, the various INPUTS and OUTPUTS are as follows:
%INPUT
%X0 is the central price
%Xmax and Xmin are the maximum and minimum of the underlying prices we want
%Delta_x is the increment in stock price
%Usin u is another way to get increment in stock prices
%OUTPUT
%Xt includes all the stock prices in the grid at each time
%rows is how many rows are in the grid
%gridS is the row in the grid where the current stock price is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the range we need so we cover all prices(Xmax and Xmin) that we want.
Universe = max(X0-Xmin, Xmax-X0)/delta_x;
N = ceil(Universe); %Making sure it is an integer
rows = 2*N + 1;%to have N on both sides of the grid
gridS = N + 1;
Xt = zeros(rows, 1);
Xt(gridS, 1) = X0;%Set the current price
%Generation of Stock Price Dynamics
for i = 1:N
    if(u == 0)
        Xt(gridS+i,1) = X0 - i*delta_x;
        Xt(gridS-i,1) = X0 + i*delta_x;
    else
        Xt(gridS+i,1) = X0*u^(-i);
        Xt(gridS-i,1) = X0*u^(i);
    end;
end;

end

