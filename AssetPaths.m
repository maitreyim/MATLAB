S0 =50;       % Price of underlying today
X = 50;       % Strike at expiry
mu = 0.04;    % expected return
sig = 0.1;    % expected vol.
r = 0.03;     % Risk free rate
dt = 1/365;   % time steps
etime = 50;   % days to expiry
T = dt*etime; % years to expiry
steps=1000;
nsims = 1000; % Number of simulated paths

% Generate potential future asset paths
% Function to generate sample paths for assets assuming geometric
% Brownian motion.
%
% S = AssetPaths(S0,mu,sig,dt,steps,nsims)
%
% Inputs: S0 - stock price
%       : mu - expected return
%       : sig - volatility
%       : dt - size of time steps
%       : steps - number of time steps to calculate
%       : nsims - number of simulation paths to generate
%
% Output: S - a matrix where each column represents a simulated
%             asset price path.
%

% calculate the drift
nu = mu - sig*sig/2;

% Generate potential paths
S = S0*[ones(1,nsims); ...
            cumprod(exp(nu*dt+sig*sqrt(dt)*randn(steps,nsims)),1)];    
        % Script to price an Asian option using a monte-carlo approach.
        % calculate the payoff for each path for a Put
PutPayoffT = max(X-mean(S),0);

% calculate the payoff for each path for a Call
CallPayoffT = max(mean(S)-X,0);

% discount back
putPrice = mean(PutPayoffT)*exp(-r*T)
callPrice = mean(CallPayoffT)*exp(-r*T) 

