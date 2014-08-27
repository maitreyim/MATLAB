tic;
% Implicit Finite Difference Method for American Put Option
% Parameter Initialization
S = 1; Smax = 100; K = 10; r = 0.10; vol = 0.30; T = 6/12;
% Grid Initialization
N = 149; % Number of Time Intervals;
M = 300; % Number of Stock Price Intervals
timestep = T/N;
stockstep = Smax/M;
% Boundary Conditions
% Payoff at Maturity
for j = 1 : M+1
f(j,N+1) = max(K - (j-1)*stockstep,0);
end;
% Stock Price = 0
for i = 1 : N+1
f(1,i) = K;
end;
% Deep Out-the-Money Put
for i = 1 : N+1
f(M+1,i) = 0;
end;
% Calculate a,b,c (parameters) for Model
for j = 1 : M-1
a(j) = 0.5*r*j*timestep - 0.5*vol^2*j^2*timestep;
b(j) = 1 + vol^2*j^2*timestep + r*timestep;
c(j) = -0.5*r*j*timestep - 0.5*vol^2*j^2*timestep;
end;
% Calculations (Solve System Ax = d)
% Set up matrix A
for j = 2 : M
A(j,j-1) = a(j-1);
A(j,j) = b(j-1);
A(j,j+1) = c(j-1);
end;
A(1,1) = 1;
A(M+1,M+1) = 1;
% Invert Matrix A
A = inv(A);
% Solve System
for i = N:-1:1
f(:,i) = A*f(:,i+1); % Calculate Value
% American Condition
for j = 1 : M+1
if f(j,i) < (K-(j-1)*stockstep)
f(j,i) = (K-(j-1)*stockstep);
end;
end;
end;
% Display Output
[call,put] = blsprice(S,K,0,T,vol,0);
disp('The Black-Scholes Price is:');
disp(put);
disp('The Finite-Difference Approximation is:');