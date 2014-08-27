clc;clear all;close all;
%Q:  
%leverage constraints is not used yet
%need to confirm the functions for Basel requirements
%%
%%%%%%%%%%%%%%%%%%%%%Assumptions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asset equity and liability weights
%assets
% Weights of cash, bond, loanbank, loancustomer
Wassets=[0.2;0.1;0.1;0.6];%sum=1
Nassets=size(Wassets,1);
RSFratio=[0,0.05,0.5,0.85];
%Weights of liability &equity
%Weights of deposit,bond,tier1,tier2
Wliab=[0.45;0.45;0.03;0.02];%sum <1
Wequity=1-sum(Wliab,1);
ASFratio=[0,1,0.85,0.7];
%Return rate of cash,bond,loanbank,loancustomer
[NUMERIC,TXT,RAW]=xlsread('C:\Users\maitreyi\Desktop\Data_set.xlsx','Data','b2:f2519');
Assetret=cell2mat(RAW);
return1=Assetret(:,1);
return2=Assetret(:,2);
return3=Assetret(:,3);
return4=Assetret(:,4);
%mean return
meanr=mean(Assetret,1)*360;
Rassets=meanr';
%covariance
for i=1:4
    for j=1:4
        assetcov=COV(Assetret(:,i),Assetret(:,j));%covariance matrix
        Covasset(i,j)=assetcov(2,2);%variance or covariance 
    end
end

% Covasset=[0.0004,0.0007,0.004,0.004;...
%           0.0007,0.002,0.003,0.003;...
%           0.004,0.003,0.006,0.006;...
%           0.004,0.003,0.006,0.006];

Eloss=[0;0.001;0.003;0.01];%based on asset, need to update%%%%%%%%%%%%%%%%%%%%
%correlation of banks
rho=0.5;
%the other bank 
Wliab2=[0.7;0.2;0.03;0.02];%the main difference between two banks
Wassets2=[0.2;0.1;0.1;0.6];
%%
%Basel III requirements
%capital adequacy ratio >= 10.5% (minimum requirement+CCB) 
%RWA calculation cash, bond,loan to bank,loan to customer
Riskweights=[0;0;0.2;1];
%RWAratio=Wequity/(Wassets'*Riskweights)>=10.5%
%equivalent to Riskweights'*Wassets<=10.5%*Wequity
A1=Riskweights';
B1=Wequity/0.105;
%leverage
%tier 1 capital is equity and disclosed reserve
leverage=(Wequity+Wliab(3,1))/1;%>=3%
%liquidity
%how about LCR ratio?
ASF=(Wequity+ASFratio*Wliab);
RSF=(RSFratio*Wassets);%NSFR=ASF/RSF>=100%
A2=RSFratio;
B2=ASF;
%not sure which is stabledeposit which is otherdeposits, P51 in Moody's slides
%%
%optimization
%CEO
UCEO= CEO(Wassets,Rassets,Wliab,Covasset,Eloss);%minus of utility function of CEO

fun = @(Wassets)CEO(Wassets,Rassets,Wliab,Covasset,Eloss);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  
A=[A1;A2];
B=[B1;B2];
Aeq=ones(1,Nassets);
Beq=1;
LB=zeros(Nassets,1);
UB=ones(Nassets,1);
[optWCEO RESNORM MaxCEO exitflag] = fmincon(fun,Wassets,A,B,Aeq,Beq,LB,UB);
MaxUCEO= -CEO(optWCEO,Rassets,Wliab,Covasset,Eloss);
%risk manager
%Umanager= manager(Wassets,Rassets,Wequity,Cov);%utility function of manager
Umanager= manager(Wassets,Rassets,Wliab,Covasset,Eloss);%the lower volatility the better

fun = @(Wassets)manager(Wassets,Rassets,Wliab,Covasset,Eloss);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  

A=[A1;A2];
B=[B1;B2];

Aeq=ones(1,Nassets);
Beq=1;
LB=zeros(Nassets,1);
UB=ones(Nassets,1);
[optWmanager RESNORM Maxmanager exitflag] = fmincon(fun,Wassets,A,B,Aeq,Beq,LB,UB);
MaxUmanager= manager(optWmanager,Rassets,Wliab,Covasset,Eloss);

%regulator
%Three policies regulator can change
RWAratio=Wequity/(Wassets'*Riskweights);%>=7%
NSFR=ASF/RSF;%>=100%
leverage=(Wequity+Wliab(3,1))/1;%>=3%
Policy=[RWAratio,NSFR,leverage];

Uregulator= regulator(Policy,Wassets,Rassets,Wliab,Wassets2,Wliab2,ASFratio,RSFratio,Covasset,Eloss,rho,Riskweights);
fun = @(Policy)regulator(Policy,Wassets,Rassets,Wliab,Wassets2,Wliab2,ASFratio,RSFratio,Covasset,Eloss,rho,Riskweights);
options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  

LB=[0,0.5,0];
UB=[0.2,5,1];

[optpolicy RESNORM Maxregulator exitflag] = fmincon(fun,Policy,[],[],[],[],LB,UB);
MinUregulator= regulator(optpolicy,Wassets,Rassets,Wliab,Wassets2,Wliab2,ASFratio,RSFratio,Covasset,Eloss,rho,Riskweights);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%CEO%%%%%
%function Utility=CEO(Wassets,Rassets,Wliab,Covasset,Eloss)
Utility=-(Wassets'*(Rassets-Eloss)-Wassets'*Covasset*Wassets);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%MANAGER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Utility=manager(Wassets,Rassets,Wliab,Covasset,Eloss)
Utility=(Wassets'*Covasset*Wassets);
end
%%%%REGULATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function utility= regulator(policy,Wassets,Rassets,Wliab,Wassets2,Wliab2,ASFratio,RSFratio,Covasset,Eloss,rho,Riskweights)
Wequity=1-sum(Wliab,1);
Wequity2=1-sum(Wliab2,1);


%Three constraints
RWAratio=policy(1);
NSFR=policy(2);
leverage=policy(3);%>=3%

%Optimization based on CEO's utility, bank 1
A1=Riskweights';
B1=RWAratio*Wequity;
%leverage
%tier 1 capital is equity and disclosed reserve
leverage=(Wequity+Wliab(3,1))/1;%>=3%
%liquidity
ASF=(Wequity+ASFratio*Wliab);
RSF=(RSFratio*Wassets);%NSFR=ASF/RSF>=100%
A2=RSFratio;
B2=ASF./NSFR;
A=[A1;A2];
B=[B1;B2];
Nassets=size(Wassets,1);
Aeq=ones(1,Nassets);
Beq=1;
LB=zeros(Nassets,1);
UB=ones(Nassets,1);
UCEO= CEO(Wassets,Rassets,Wliab,Covasset,Eloss);%minus of utility function of CEO
fun = @(Wassets)CEO(Wassets,Rassets,Wliab,Covasset,Eloss);
options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  
[optW RESNORM Maxregulator exitflag] = fmincon(fun,Wassets,A,B,Aeq,Beq,LB,UB);
%Optimization based on CEO's utility, bank 2
A3=Riskweights';
B3=RWAratio*Wequity2;
%leverage
%tier 1 capital is equity and disclosed reserve
leverage2=(Wequity2+Wliab2(3,1))/1;%>=3%
%liquidity
ASF2=(Wequity2+ASFratio*Wliab2);
RSF2=(RSFratio*Wassets2);%NSFR=ASF/RSF>=100%
A4=RSFratio;
B4=ASF2./NSFR;
A=[A3;A4];
B=[B3;B4];
Nassets2=size(Wassets2,1);
Aeq2=ones(1,Nassets2);
Beq2=1;
LB2=zeros(Nassets2,1);
UB2=ones(Nassets2,1);
UCEO2= CEO(Wassets2,Rassets,Wliab2,Covasset,Eloss);%minus of utility function of CEO
fun = @(Wassets2)CEO(Wassets2,Rassets,Wliab2,Covasset,Eloss);
options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  
[optW2 RESNORM Maxregulator exitflag] = fmincon(fun,Wassets2,A,B,Aeq2,Beq2,LB2,UB2);
%minimize rho*sigmaA*sigmaB
utility=(optW'*Covasset*optW)*(optW2'*Covasset*optW2);
optW
optW2
end%
