clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%Assumptions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asset equity and liability weights
%assets
% Weights of cash, bond, retailloan, commloan, bankloan
Wassets=[20;480;1000;2500;500];%amount
Totalassets=sum(Wassets,1);
RSFratio=[0,0.05,0.85,0.85,0.85];
%Weights of liability &equity
Wliab=[3000;200;800];%amount
Wequity=sum(Wassets,1)-sum(Wliab,1);
ASFratio = [85,70,85];

%Return rate of assets
%mean return
Rassets=[0;0.01;0.07;0.05;0.03];
Rliab = [0.03;0.04;0.05];
%covariance
 Covasset=[0,0,0,0,0;...
           0,0.01^2,0,0,0;...
           0,0,0.09^2,0,0;...
           0,0,0,0.06^2,0;...
           0,0,0,0,0.02^2];
       
%   Covasset=[0,0,0,0,0;...
%            0,0.00,0,0,0;...
%            0,0,0.00,0,0;...
%            0,0,0,0.00,0;...
%            0,0,0,0,0.00];

Eloss=[0;0;0.05;0.09;0.01];%based on asset
%correlation of banks
rho=0.5;
%the other bank 
Wliab2=[3000;200;800];
Wassets2=[200;300;1000;2500;500];
Wequity2=sum(Wassets2,1)-sum(Wliab2,1);
%%
%Basel III requirements, assume equity is fixed
%capital adequacy ratio >= 10.5% (minimum requirement+CCB) 
%RWA calculation cash, bond,loan to bank,loan to customer
Riskweights=[0;0;1;1;0.2;0;0;0];%weight for asset and liability(0)
%RWAratio=Wequity/(Wassets'*Riskweights)>=10.5%
%equivalent to Riskweights'*Wassets<=10.5%*Wequity
A1=Riskweights';
B1=(Wequity/0.06);
%leverage
%tier 1 capital is equity 
%leverage=(Wequity)/Totalassets>=3%
A2=[1,1,1,1,1,0,0,0];
B2= (Wequity/0.03);
%liquidity
ASF=(Wequity+ASFratio*Wliab);
RSF=(RSFratio*Wassets);%NSFR=ASF/RSF>=100%
A3=[RSFratio,-1.*ASFratio];
B3=Wequity;
%%
% %optimization
% %CEO
X=[Wassets;Wliab];
X2=[Wassets2;Wliab2];
%UCEO= CEO(X,Wequity,Rassets,Rliab);%minus of utility function of CEO

fun = @(X)CEO(X,Wequity,Rassets,Rliab,Covasset);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  

A=[A1;A2;A3];
B=[B1;B2;B3];

Aeq=[1,1,1,1,1,0,0,0;
    0,0,0,0,0,1,1,1];
Beq=[Totalassets;Totalassets-Wequity];
LB=zeros(8,1);
UB=10000.*ones(8,1);
[optWCEO RESNORM MaxCEO exitflag] = fmincon(fun,X,A,B,Aeq,Beq,LB,UB);
%MaxUCEO= -CEO(optWCEO,Wequity,Rassets,Covasset,Eloss);
% %risk manager
% 
%Umanager= manager(Wassets,Rassets,Wequity,Cov);%utility function of manager
%Umanager= manager(X,Wequity,Rassets,Rliab,Covasset);%the lower volatility the better

fun = @(X) manager(X,Wequity,Rassets,Rliab,Covasset);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  

A=[A1;A2;A3];
B=[B1;B2;B3];


Aeq=[1,1,1,1,1,0,0,0;
    0,0,0,0,0,1,1,1];
Beq=[Totalassets;Totalassets-Wequity];
LB=zeros(8,1);
UB=10000.*ones(8,1);
[optWmanager RESNORM Maxmanager exitflag] = fmincon(fun,X,A,B,Aeq,Beq,LB,UB);
MaxUmanager= manager(optWmanager,Wequity,Rassets,Covasset,Eloss);
