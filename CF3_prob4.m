%Problem 4: Pricing European Put Option by Monte Carlo Simulation
% European Call Parameters As Given with S0,X (Exercise Price),r(interest rate),
%T,sigma(volatility)        
N=5000; % # of obs
r=0.04;% interest rate as given in the question
T=0.5;
sigma=0.25;% volatility as given in the question
S0=20;% initial stock price as given in the question
X=20;% Exercise Price as given in the question
W =sqrt(T)*randn(N,1); %generation of standard wiener process
% Stock Price Dynamics with Geometric Brwonian Motion
Sp=S0*exp(sigma.*W+(r-(sigma^2)/2)*T);% Creating stock price dynamics
%European Put Price with W
dp=exp(-r*T)*max(X-Sp,0); 
c1=mean(dp);% European put price
%elseif dp=0;
%end
 Probability1= mean( (X-Sp) > 0 );% computing probability that European put price will exercise in the money
% % Antithetic Variates Technique
Sm=S0*exp(sigma.*W*(-1)+(r-(sigma^2)/2)*T);% Creating stock price dynamics
% European Put Price with W
dm=exp(-r*T)*max(X-Sm,0); 
c2=mean(dm);% European put price
Probability2= mean( (X-Sm) > 0 );% computing probability that European put price will exercise in the money



