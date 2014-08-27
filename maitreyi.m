% Question 4: Comparison of European Vs. American Put Pricing
% Main Script to call functions created for European & American Put
clear all

r = .03;% Risk free interest rate as given in the question
T = 1;% Total Time
sigma = .3; % Volatility
K = 100; % Exercise or Strike Rate

S0_set = (80:4:120)'; % Set of Stock Price to vary from $80 to $120 with an increment of $4

APs = nan(size(S0_set, 1), 1); % Preallocating memory space for American Puts
EPs = nan(size(S0_set, 1), 1);%Preallocating memory space for European Puts

for i = 1:size(S0_set, 1) % Looping over stock price
    EPs(i) = CF4_EP(S0_set(i),T,r,sigma,K );% European Put Prices
    APs(i) = CF4_AP(S0_set(i),T,r,sigma,K );% American Put Prices
end
% Plotting the European Vs American
plot(S0_set, [EPs, APs])
Legend('European Put', 'American Put');
Xlabel('Stock Price Dynamics');
Ylabel('Put Option Price');
Title('European Vs American Put Price comparison Using Binomial');


