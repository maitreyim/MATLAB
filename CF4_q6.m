%Question 6: Using Halton's Low-Discrepancy Sequences to price European call options 
%Halton's Low-Discrepancy Sequence
clear all;
N_vec = 100; % Length of vector
seed = 98242;% Use the seed
trials = 2000;% # of iterations
haltonSeqB = zeros(N_vec/2 , 2);% Create The Vector
b1=2; % First Halton Series of Base 1
b2=7; % First Halton Series of Base 1
r = 0.05; % interest rate as given
sigma = 0.24; % volatility rate as given
S0=32;% Initial Stock Price as given
K=30; % Strike or Exercise Rate as given
T=0.5; % Time as given in the problem

seqBase_b1 = GetHalton(N_vec/2 , b1); % Use the function GetHalton (base b1) as defined
seqBase_b2 = GetHalton(N_vec/2 , b2);% Use the function GetHalton (base b2) as defined

for i = 1:N_vec/2 % looping Over the length of the vector
    haltonSeqB(i, 1) = seqBase_b1(i);
    haltonSeqB(i, 2) = seqBase_b2(i);
end; 
%Box-Muller to generate Normal Distribution
  z1(:,1) = sqrt(-2*log(haltonSeqB(:, 1))).*cos(2*pi*haltonSeqB(:, 2));
  z2(:,1)= sqrt(-2*log(haltonSeqB(:, 1))).*sin(2*pi*haltonSeqB(:, 2));
z=([z1; z2]);% Stacking the Normals
% Call Pricing
call= (exp(-r*T))*max(((S0*(exp((r-((sigma^2)/2)*T)+sigma*z*sqrt(T))))-K),0);
% Final Call Value
mean(call);