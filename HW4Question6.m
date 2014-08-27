function HW4Question6()
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
T = 0.5; %Time to maturity
r = 0.05;  %Interest rate
sigma = 0.24; %Volatility
S0 = 32; %Stock price
K = 30; %Strike price
%Vector with different amount of steps in the tree
nVec = [10, 15, 20, 40, 70, 80, 100, 200, 500];

%a)
nOpt = zeros(1, 9);
nReal = zeros(1, 9);
for i = 1:9
    n = nVec(i); %Amount of steps in the tree
    delta = T / n; %Time interval
    d = exp(-sigma*sqrt(3*delta));
    u = 1/d;
    pd = (r*delta*(1-u)+(r*delta)^2+sigma^2*delta)/((u-d)*(1-d));
    pu = (r*delta*(1-d)+(r*delta)^2+sigma^2*delta)/((u-d)*(u-1));
    %Get the price of a call option using Trinomial tree
    nOpt(i) = PriceAnOptionTrinomial('E', 'C', S0, K, r, sigma, T, n, u, d, pu, pd, 'R');
    %Get the price of a call option using Black-Scholes
    nReal(i) = blsprice(S0, K, r, T, sigma);
end;

%b)
nOptLog = zeros(1, 9);
for i = 1:9
    n = nVec(i); %Amount of steps in the tree
    delta = T / n; %Time interval
    d = -sigma*sqrt(3*delta);
    u = sigma*sqrt(3*delta);
    pd = 0.5*((sigma^2*delta+(r-sigma^2/2)^2*delta^2)/u^2-(r-sigma^2/2)*delta/u);
    pu = 0.5*((sigma^2*delta+(r-sigma^2/2)^2*delta^2)/u^2+(r-sigma^2/2)*delta/u);
    %Get the price of a call option using Trinomial tree and log prices
    nOptLog(i) = PriceAnOptionTrinomial('E', 'C', log(S0), log(K), r, sigma, T, n, u, d, pu, pd, 'L');
end;
%Plot trinomial compared to Black-Scholes
figure(1)
plot(nVec, nOpt, '-r', nVec, nReal, '-k', 'linewidth',2);
legend(3,'Trinomial', 'Black-Scholes');
%Plot Log
figure(2)
plot(nVec, nOptLog, '-r', 'linewidth',2);
legend(3,'Trinomial-Log');
end

