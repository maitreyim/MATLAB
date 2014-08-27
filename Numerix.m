function [ Price ] = Numerix(r0, sigma, kappa, aveR, T, n)
%Function Numerix 
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
R=0.08;
Pv0=100000;
Month_fac=repmat([0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98],1,30);
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
PV = zeros(n,1);
PV(1) = Pv0;
TPP = zeros(n,1);
r=R/12;
for j=1:1000
    Z = randn(n, 1);%Create random numbers for each path
    for i = 2:n%To create each point in the path
        deltaR = kappa*(aveR - rt(i-1))*deltaT + sigma*sqrt(rt(i-1))*sqrt(deltaT)*Z(i-1,1);
        rt(i) = rt(i-1) + deltaR;
        RI=0.28 + 0.14*atan(-8.57+430*(R-rt(i-1)));%Interest rate Factor
        BU=0.3+0.7*PV(i-1)/Pv0;%BurnOut Factor
        SG=min(1,i/30); % Seasonality Factor
        SY=Month_fac(i);% Seasoning Factor
        CPR=RI*BU*SG*SY;% Conditional Prepayment Modeling
        term = 1/( 1 - (1+r)^(-n + (i-1))) - 1;
        TPP(i) = PV(i-1) * ( r*term + ( 1 - r*term) * (1 - (1 - CPR)^(1/12)));
        c(i) =  PV(i-1) * ( r*(term +1) + ( 1 - r*term) * (1 - (1 - CPR)^(1/12)));
        PV(i)=PV(i-1)-TPP(i);
            end;

end
   discount = (1./(cumprod(1+rt))).^(1/12);
   Price=sum(c.*discount')% MBS Price
    
    
    
    
    
    
    
    
    

