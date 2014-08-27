function [ curve ] = OneFactorCIR( r0, sigma, kappa, aveR, T, n, M )
%Function to generate stochastic interest rate using CIR model
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
r = zeros(1,M);%Average rate of each path

for j = 1:M%To create each path
    Z = randn(n, 1);%Create random numbers for each path
    for i = 1:n-1%To create each point in the path
        deltaR = kappa*(aveR - rt(i))*deltaT + sigma*sqrt(rt(i))*sqrt(deltaT)*Z(i,1);
        rt(i+1) = rt(i) + deltaR;
    end;
    if(M ~= 1)%If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;
curve = r;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ curve ] = OneFactorVasicek( r0, sigma, kappa, aveR, T, n, M )
%Function to generate stochastic interest rate using Vasicek model
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
r = zeros(1,M);%Average rate of each path

for j = 1:M %To create each path
    Z = randn(n,1); %Create random numbers for each path
    for i = 1:n-1 %To create each point in the path
        deltaR = kappa*(aveR - rt(i))*deltaT + sigma*sqrt(deltaT)*Z(i,1);
        rt(i+1) = rt(i) + deltaR;
    end;
    if(M ~= 1) %If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;

curve = r;
end