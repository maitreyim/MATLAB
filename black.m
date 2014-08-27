%************* Black's Approximation ********************
% For American Call Options 
% clear all; clc; close all;
m=[1.0 4.0 7.0]; T=m/12; r =0.05; si = 0.3; K = [35 40 45]; %define stock process variables

S0 = 40; % define the initial Stock prices
D={'Cash' 'tau1' 'Frequency'  'Quarterly-Yield'   'CashOrYield' '#BinomialSteps' 'n'      'T';
    0.5     0.5   3               0.0125          'N'             1              0      T(1)};
ctr=0;
for i=1:3 % strike price loop
    
    for j=1:3 % time to maturity loop
        ctr=ctr+1;
        [ option ] = BlackApprox( r, si, S0, K(i), T(j), D);
        Opt(ctr,1)=K(i);Opt(ctr,2)=T(j);Opt(ctr,3)=option;
    end
end
