function CF9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTATIONAL FINANCE PROJECT #9
%MORTGAGE BACKED SECURITIES MODELING%
%%%SUBMITTED BY MAITREYI MANDAL%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%Question 1A%%%%%%%%%%
%%PARAMETER INITIALIZATION%%%%%%
% aveR=0.08;%Rbar as given
% sigma=0.09;%Volatility as given
% R=0.08; 
% r0=0.078;
% T=30;%# of years
% n=360;% no. of months
% kappa_set=0.3:0.1:0.9;% set of kappa need to be used
% for i = 1:size(kappa_set,2)
% Price(i)=Numerix( r0, sigma, kappa_set(i), aveR, T, n);
% end
% subplot(2,2,1)
% % Plotting the graph w.r.t. kappa set
% plot(kappa_set,Price)
% legend('Price','Kappa');
% Title('MBS Vs Kappa');
% Xlabel('Kappa');
% Ylabel('MBS Price(Numerix)')
% % % % %%%%%QUESTION 1B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma=0.09;
% kappa=0.6;
% R=0.08;
% T=30;
% n=360;
% r0=0.078;
% aveR_set=0.03:0.01:0.09;
% for i = 1:size(aveR_set,2)
%      
% Price(i)=Numerix( r0, sigma, kappa, aveR_set(i), T, n);
% end
% subplot(2,2,2)
% Plotting the graph w.r.t. Rbar set
% plot(aveR_set,Price)
% legend('Price','Rbar');
% Title('MBS Vs Rbar(Numerix)');
% Xlabel('Rbar');
% Ylabel('MBS Price(Numerix)')
% % % %%%%%%%%%%%QUESTION 2A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveR=0.08;
sigma=0.09;
R=0.08;
T=30;
n=360;
r0=0.078;
kappa_set=0.3:0.1:0.9;
for i = 1:size(kappa_set,2)
      
Price(i)=PSA( r0, sigma, kappa_set(i), aveR, T, n);
    
end
% Plotting the graph w.r.t. kappa set
subplot(2,2,3)
plot(kappa_set,Price)
legend('Price','Kappa');
Title('MBS Vs Kappa (PSA)');
Xlabel('Kappa');
Ylabel('MBS Price(PSA)')
% % % % %%%%%%%%%%%QUESTION 2B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r0=0.078;
% sigma=0.09;
% R=0.08;
% T=30;
% n=360;
% aveR_set=0.03:0.01:0.09;
% for i = 1:size(aveR_set,2)
%       
%    Price(i)=PSA( r0,sigma, kappa, aveR_set(i), T, n);
%    
% end
% % Plotting the graph w.r.t. Rbar set
% subplot(2,2,4)
% plot(aveR_set,Price)
% legend('Price','Rbar');
% Title('MBS Vs Rbar(PSA)');
% Xlabel('Rbar');
% Ylabel('MBS Price(PSA)')
end
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
M = 10000;
for j=1:M
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
end
   function [ Price ] = PSA( r0, sigma, kappa, aveR, T, n)
%Function to generate stochastic interest rate using CIR model
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
R=0.08;
Pv0=100000;
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
PV = zeros(n,1);
PV(1) = Pv0;
TPP = zeros(n,1);
r=R/12;
c = zeros(n,1);
for j = 1:1000%To create each path
    Z = randn(n, 1);%Create random numbers for each path
    for i = 2:n-1%To create each point in the path
        deltaR = kappa*(aveR - rt(i-1))*deltaT + sigma*sqrt(rt(i-1))*sqrt(deltaT)*Z(i-1,1);
        rt(i) = rt(i-1) + deltaR;
        term = 1/( 1 - (1+r)^(-n + (i-1))) - 1;
        CPR=min(0.06,(i-1)*0.002);
        TPP(i) = PV(i-1) * ( r*term + ( 1 - r*term) * (1 - (1 - CPR)^(1/12)));
        c(i) =  PV(i-1) * ( r*(term +1) + ( 1 - r*term) * (1 - (1 - CPR)^(1/12)));
        PV(i)=PV(i-1)-TPP(i);
        
    end;
%     figure(2)
%     plotyy(1:n, TPP, 1:n, rt)
end
    
    discount = (1./(cumprod(1+rt))).^(1/12);
    Price=sum( c .* discount )% MBS Price Using PSA
   end
