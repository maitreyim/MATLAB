%Problem 4: Pricing European Call Option by Monte Carlo Simulation
% European Call Parameters As Given with S0,X (Exercise Price),r(interest rate),
%T,sigma(volatility)        
N=5000000;
T=10;
r=0.04;
sigma=0.2;
S0=88;
X=100;
W =sqrt(T)*randn(N,1);
% Stock Price Dynamics with Geometric Brwonian Motion
Sp=S0*exp(sigma.*W+(r-(sigma^2)/2)*T);
%European Call Price with W
dp=exp(-r*T)*max(Sp-X,0);
c1=mean(dp);% European Call Price given by c1 using W()
% Using Antithetic Variates of Variance Reduction Method
Sm=S0*exp(sigma.*W*(-1)+(r-(sigma^2)/2)*T);% Creating Sminus with W(-1)
%European Call Price with W(-1)
dm=exp(-r*T)*max(Sm-X,0);
d = (dp + dm)/2;
c2=mean(d);% European Call Price given by 
%Antithetic Variates of Variance reduction Technique
%European Call Price Using built-in MATLAB command of blsprice
[call]=blsprice(88,100,0.04,10,0.2);










