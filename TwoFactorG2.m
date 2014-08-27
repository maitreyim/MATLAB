function [ curve ] = TwoFactorG2( x0, y0, r0, phi0, rho, sigma, nu, a, b, T, n, M )
%Function to generate stochastic interest rate using G2++ model
%This is a two factor model where the interest rates are generated from two
%factors.
%x0 and y0 are the initial values of the two factors
%phiO is the initial value of the shift funtion phit.
%r0 is the initial interest rate, sigma and nu are the volatilities for xt and yt,
%a and b are the drifts for xt and yt
%rho is the correlation between the two factors
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths. n is the number of steps.

deltaT = T/n;%Size of each time step
%Vectors to capture all the points in one path of interest rate
rt = zeros(n,1); xt = zeros(n,1); yt = zeros(n,1); phit = zeros(n,1);
%The first point in the path
rt(1) = r0; xt(1) = x0; yt(1) = y0; phit(1) = phi0;
r = zeros(1,M);%Average rate of each path

for j = 1:M%To create each path
    Z = randn(n,2);%Create random numbers for each path, both xt and yt
    for i = 1:n-1%To create each point in the path
        %We need to find the changes in xt, yt and phit
        dx = -a*xt(i)*deltaT+sigma*sqrt(deltaT)*Z(i,1);
        dy = -b*yt(i)*deltaT+nu*sqrt(deltaT)*(rho*Z(i,1)+sqrt(1-rho^2)*Z(i,2));
        dphi = 0;
        %Get the next point for xt, yt and phit
        xt(i+1) = xt(i) + dx;
        yt(i+1) = yt(i) + dy;
        phit(i+1) = phit(i) + dphi;%phi is a shift function
        %Get each interest rate point, r, from xt, yt and phit
        rt(i+1) = xt(i+1) + yt(i+1) + phit(i+1);
    end;
    if(M ~= 1)%If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;

curve = r;
end

