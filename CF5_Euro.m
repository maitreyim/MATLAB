%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quantitative Finance  Assignment #5
%Submitted By: Maitreyi Mandal
%Question 3 a and 3 b
%European and American Forward Start Put Pricing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 3a: Forward Start European Put Pricing
S0=65;%Initial Stcok Price as asked in the question
sigma=.2;% Volatility as provided
r=.06; % Risk Free Interest Rate
t=.2; %small t
T=1;% Big T, entire Time Interval
N=1000000; %# of simulations
St=S0*exp((r-(sigma^2)/2)*t+sqrt(t)*sigma*randn(N,1));% Stock Price at St

% % [Call,Put] = blsprice(Price, Strike, Rate, Time, Volatility)
%  
 [C,P] = blsprice(St,St,r,T-t,sigma);% Just using in-built 
 %blsprice to compute Forward Start European Put Option Price 
 mean(P) % Forward Start European Put Price
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now Comes the Harder Part: Forward Start American Put Price
Smin=min(St);% Checking what the minimum stock price is
Smax=max(St);% Checking what the maximum stock price is
nbin=200;
Amputstrike=linspace(Smin,Smax,nbin)'; % Creating a grid with nbin for Strike Prices
Amputprice=nan(nbin,1); % Defining the American Put Price
for i=1:nbin
    Amputprice(i,1)=CF4_AP(Amputstrike(i,1),T-t,r,sigma,Amputstrike(i,1) );
end
% CF4_AP was previously written for assignmnet 4 
counts=nan(nbin-1,1);
for j = 2:nbin
    counts(j-1,1)=sum((St>=Amputstrike(j-1,1))&((St<=Amputstrike(j,1))));
end
% summing up how many times St>= and <= Amputstrike
binmidpoints = (Amputprice(1:(end-1),1) + Amputprice(2:end, 1))/2;

weights = counts / sum(counts); % Creating weights
sum( weights .* binmidpoints )% Forward Start American Put Price


%plot(Amputstrike, weights)



% 
