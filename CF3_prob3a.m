%Problem 3, Part a
%European Call Function Approximation
function [E_call,E_call2] = CF3_prob3a(S0,T, X, si, r)
if si<0
    error( 'Volatility cannot be negative')
end
W =sqrt(T)*randn(1000,1);
% Stock Price Dynamics with Geometric Brwonian Motion
Sp=S0*exp(si.*W+(r-(si^2)/2)*T);
%European Call Price with W
dp=exp(-r*T)*max(Sp-X,0);
E_call= mean(dp);% European Call Price given by c1 using W()
% Using Antithetic Variates of Variance Reduction Method
Sm=S0*exp(si.*W*(-1)+(r-(si^2)/2)*T);% Creating Sminus with W(-1)
%European Call Price with W(-1)
dm=exp(-r*T)*max(Sm-X,0);
d = (dp + dm)/2;
E_call2=mean(d);% European Call Price given by 
% Antithetic Variates of Variance reduction Technique
end

