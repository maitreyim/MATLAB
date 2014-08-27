%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computational Finance Homework#5
%Question#1A: Estimation of American Put Option by Least-Square 
%MonteCarlo Using Different Basis Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Parameters Given in the question
T = 0.5;% Time Period
r = .06;% Interest Rate
S0 = 36;% Initial Stock Price As Provided
K = 40; % Strike Price As Provided
sigma = .2;% Volatility As Provided
N = 200; 
M = 50000;% # of paths
z1=randn(M,N);% Gdenerating normal random number
z2=randn(M,N)*(-1);% Antithetic Variable

h = T / N;% Discretization of time step

diffusionp = [zeros(M, 1) sigma * sqrt(h) * cumsum( z1, 2)];% Calculating the diffusion component
diffusionm=[zeros(M, 1) sigma * sqrt(h) * cumsum( z2, 2)];%antithetic var
drift     = [zeros(M, 1) repmat((r-sigma^2/2)*h*(1:N), M, 1) ];% Calculating the drift Component
% Capturing Stock Price Dynamics
Sp = S0*exp(drift + diffusionp);
Sm = S0*exp(drift + diffusionm);
S=[Sp;Sm];
dmat = repmat( exp(-r*h*(0:N)), M*2, 1);
ind = zeros(M*2, N+1);% Index Generation
ind(:, end) = ( (K - S(:, end))>0 );
val = max(0, K - S);
% Running the Itearation
n = N;
for n = N:(-1):1
  if( mod(n, 50) == 0 )
    disp(['Iteration ' num2str(n) ' complete'])% Just a check to 
    %see if the code is running
  end
  V = val(:,n);
  itm = (V > 0);
  % Continuation Value
  
  CV = sum( val .* dmat .* ind, 2);
  
  bases = Hermite(S(:,n), 3);
  
  Y = CV(itm,:);
  X = bases(itm,:);
  beta = regress(Y,X);% Running the regression
  ECV = bases * beta;
  
  ex = (itm & ( ECV < V) );
  ind(ex,:) =0;
  ind(ex, n) = 1;

end

mean(  sum( val .* dmat .* ind, 2) ) % Finally Calculating the American Put Pricing











