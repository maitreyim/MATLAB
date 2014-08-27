%Question 2
clear all
S0    = 615;% Current GOOGLE Price taken from Yahoo! Finance
T     = 9/12;% Time Period
r     = .02; % Risk-Free rate
sigma = .38;% Volatility as obtained from Historical 60 months of Google Stock Price
K=round(1.1*S0); % Setting the Strike/Exercise Price
M = 200; % No of Steps taken
delta = T / M; % each time step

% Calculating Google Call Price Using 1b method
u=exp(r*delta)*(1+sqrt(exp((sigma^2)*delta)-1)); 
d=exp(r*delta)*(1-sqrt(exp((sigma^2)*delta)-1));
p=0.5;

Goog_C = max(S0 * ( (d.^(0:M)) .* (u.^(M:(-1):0)) )' - K, 0);
for m = 1:M
    Goog_C = 1/(1 + r * delta ) * ( Goog_C(1:(end-1)) * p + Goog_C(2:end) * (1-p) );
end

Goog_C












