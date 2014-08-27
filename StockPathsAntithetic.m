function S = StockPathsAntithetic(S0,r,si,dt,steps,nsims)
% Function to generate sample paths using antithetic variates
% When rnd is a vector of random numbers used to generate
% one simulated path then -rnd is used to generate another
% path.
%
% S = StockPathsAntithetic(S0,r,si,dt,steps,nsims)
%
% Inputs: S0 - stock price
%       : mu - expected return
%       : si - volatility
%       : dt - size of time steps
%       : steps - number of time steps to calculate
%       : nsims - (even) number of simulation paths to generate
%
% Output: S - a matrix where each column represents a simulated
%             asset price path.
%
% Notes: This code focuses on details of the implementation of the
%        Monte-Carlo algorithm.

% calculate the drift
dr = r - si*si/2;

% generate random samples
rnd = randn(steps,nsims/2);

% Generate potential paths
S = S0*[ones(1,nsims); ...
            cumprod(exp(dr*dt+si*sqrt(dt)*[rnd -rnd]),1)];

