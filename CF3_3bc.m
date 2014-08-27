%Part B: Write a code to compute 
% set the parameters
S0    = 100;
T     = 1;
r     = .04;
sigma = .25;
K     = 95;
d1 = (log(S0/K) + (r+sigma^2/2)*T)/(sigma*sqrt(T));%d1 defined
d2 = (log(S0/K) + (r-sigma^2/2)*T)/(sigma*sqrt(T));%d2 defined
Nd1=getN(d1);%Nd1 defined
Nd2=getN(d2);%Nd2 defined
% Get the European Call Price
C= S0 * Nd1- exp(-r*T) * K * Nd2;%Call price
%Question 3, Part c)
%OPTION GREEKS
call = blsprice(90, 105, 0.04, 5, 0.2); 
increment = 0.0005; % step size as asked
r = 0.04; % interest rate
sigma = 0.25; % volatility
X = 20; % initial stock price
T = 0.5; % time
delta = zeros(1, 11);
gamma = zeros(1, 11);
theta = zeros(1, 11);
vega = zeros(1, 11);
rho = zeros(1, 11);

for i = 15:1:25
    S0 = i;
    %Computing Delta: the First Option Greek, which measures sensitivity of
    %change in call price w.r.t change in stock price, this is equal to
    %N(d1)
    S1 = S0;
    S2 = S0 + increment;
    First_S1 = blsprice(S1, X, r, T, sigma);
    Next_S2 = blsprice(S2, X, r, T, sigma); 
    delta1 = (Next_S2 - First_S1) / increment;
    delta(i - 14) = delta1;
    %Using our Defined Function that we created called getN()
    d1 = (log(S0/X) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T));
	d2 = d1 - sigma * sqrt(T);
    N1 = getN(d1);
    N2 = getN(d2);
    %Computing Gamma: The Second Greek, which measures change of delta with
    %respect to change in stock price, this is actually a pdf
    S3 = S0 + 2 * increment;
    Last_S3 = blsprice(S3, X, r, T, sigma);
    delta2 = ( Last_S3  -  Next_S2) / increment;
    gamma(i - 14) = (delta2 - delta1) / increment;
       
    %computing Theta: the third option greek, which measures sensitivity of
    %cal;l price w.r.t change in time
    T1 = T;
    T2 = T + increment;
    First_T1 = blsprice(S0, X, r, T1, sigma);
    Next_T2 = blsprice(S0, X, r, T2, sigma); 
    theta(i - 14) = ( Next_T2 - First_T1) / increment;
    %Computing Vega: the fourth greek, which measures sensitivity of call
    %price w.r.t. change in volatility
    sigma1 = sigma;
    sigma2 = sigma + increment;
    First_Sigma1 = blsprice(S0, X, r, T, sigma1);
    Next_Sigma2 = blsprice(S0, X, r, T, sigma2); 
    vega(i - 14) = (Next_Sigma2 - First_Sigma1) / increment;
     
    %Computing Rho: the fifth greek, which measures sensitivity of call
    %price w.r.t. change in interest rate
    r1 = r;
    r2 = r + increment;
    First_R1 = blsprice(S0, X, r1, T, sigma);
    Next_R2 = blsprice(S0, X, r2, T, sigma); 
    rho(i - 14) = ( Next_R2- First_R1) /increment;
       end;

plot(15:1:25, delta, '-b', 15:1:25, gamma, '-g', 15:1:25, theta, '-r', 15:1:25, vega, '-c', 15:1:25, rho, '-k', 'linewidth',3);
legend('Delta', 'Gamma', 'Theta', 'Vega', 'Rho');
Title('Option Greeks');
Xlabel('StockPrice Dynamics')
Ylabel('Greek Dynamics')



