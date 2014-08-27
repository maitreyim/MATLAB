function PathAvgOption
% %PathAvgOption Submission for Term Paper Initial function to replicate
% %   Table 1 results from Kemna-Vorst paper of 1989
% % Average path options - Asian Options
% %% Define Headers of Table
Price{1,1}='TABLE 1';
Price{2,1}='Monte Carlo simulation results';
ctr=3;
Price{ctr,1}={'sigma'}; Price{ctr,2}={'K'};
Price{ctr,3}={'Ca'}; Price{ctr,4}={'C(rnd)b'}; Price{ctr,6}={'C(rnd) red. var c'};
% %% Define input parameters
r=[.03 .05 .07]; si=[.2 .3 .4]; K=[35 40 45]; T=4/12; S0=40;%Cox and Rubinstein 1985 paper p. 216 values
% %% Define Monte Carlo parameters
N=10000;%number of simulations
n=T*12*22-1;%number of time steps based on trading days in T*12 months
% %% build the table
for i=1:3 % rate loop
    ctr=ctr+1;
    Price{ctr,1}={'r='};Price{ctr,2}= 1+r(i);
    for j=1:3 % volatility/sigma loop
        for k=1:3 % Strike Price Loop
            ctr=ctr+1;
            Price{ctr,1}=si(j);
            Price{ctr,2}=K(k);
            [Ca Cb siCb CcRedVar siCc]=getPathAvgOptPrc(N,n,T,S0,r(i),si(j),K(k));
            Price{ctr,3}=Ca; 
            Price{ctr,4}=Cb; Price{ctr,5}= siCb;
            Price{ctr,6}=CcRedVar; Price{ctr,7}=siCc;
        end
    end
end
 Price %print Table 1 results
end

%% function to get Path Avergae Option price
function [Ca Cb siCb CcRedVar siCc]=getPathAvgOptPrc(N,n,T,S0,r,si,K)
% Generate potential future asset paths
% Function to generate sample paths for assets assuming geometric
% Brownian motion.
%N      Number of simulated paths
%n      number of time steps
%S0     Price of underlying today
%K      Strike at expiry
%si     expected vol.
%r      Risk free rate
%T      Time to expiry in years

nu = r - si*si/2;% calculate the drift
dt = T/n;
% Generate potential paths
S = S0*[ones(1,N); ...
            cumprod(exp(nu*dt+si*sqrt(dt)*randn(n,N)),1)];    
        % Script to price an Asian option using a monte-carlo approach.

% Generate potential paths for log stock price
% This code is not used in the submission as I am using the normal Stock
% price generation process which gives same results as log Stock process. 
logS = log(S0)*ones(1,N);
logSt = bsxfun(@times,logS, cumsum(nu*dt + si*sqrt(dt)*rand(n,N),1));    
Sexp = exp(logSt);
        % Script to price an Asian option using a monte-carlo approach.
        

        
% calculate the payoff for each path for a standard call option
stdCall = max(S(end,:)-K,0);

% calculate the payoff for each path for an asian/artihmetic average Call
avgCall = max(mean(S)-K,0);

% calculate the payoff for each path for a geomteric average Call
geoCall = max(geomean(S)-K,0);

%calculate geometric average option using the analytic solution
nuMinus=r - si*si/6; nuPlus=r + si*si/6;
dstar=0.5*nuMinus*T;
d=(log(S0/K)+0.5*(nuPlus)*T)/(si*sqrt(T/3));
geoCallAnly=exp(dstar)*S0*normcdf(d,0,1) - K*normcdf(d-si*sqrt(T/3));
avgCallRedVar=avgCall-geoCall+geoCallAnly;

% discount back
CaMC = mean(stdCall)*exp(-r*T); % not used in code. I use Black-Scholes instead
Cb = mean(avgCall)*exp(-r*T);
siCb = STD(avgCall)/sqrt(N);

CcRedVar = mean(avgCallRedVar)*exp(-r*T);
siCc = STD(avgCallRedVar)/sqrt(N);

Ca=blsprice(S0,K,r,T,si); % BS price for standard call stock option
end