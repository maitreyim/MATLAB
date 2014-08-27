clc;clear all;close all;
%Q:  
%leverage constraints is not used yet
%need to confirm the functions for Basel requirements
%%
%%%%%%%%%%%%%%%%%%%%%Assumptions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asset equity and liability weights
%assets
% Weights of cash, bond, loanbank, loancustomer
Wassets=[0.2;0.1;0.1;0.5;0.1];%sum=1
Nassets=size(Wassets,1);
RSFratio=[0,0.05,0.5,0.85,0.9];
%Weights of liability &equity
%Weights of deposit,bond,tier1,tier2
Wliab=[0.45;0.45;0.03;0.02;0.04];%sum <1
Wequity=1-sum(Wliab,1);
ASFratio=[0,1,0.85,0.7,0.5];
%Return rate of cash,bond,loanbank,loancustomer
[NUMERIC,TXT,RAW]=xlsread('C:\Users\maitreyi\Desktop\Data_set.xlsx','Data','b2:f2518');
Assetret=cell2mat(RAW);
return1=Assetret(:,1);
return2=Assetret(:,2);
return3=Assetret(:,3);
return4=Assetret(:,4);
return5=Assetret(:,5);
%mean return
meanr=mean(Assetret,1)*360;
Rassets=meanr';
V=cov(Assetret);
Eloss=[0;0.001;0.003;0.01;0.01];%based on asset, need to update%%%%%%%%%%%%%%%%%%%%
%correlation of banks
rho=0.5;
%the other bank 
Wliab2=[0.7;0.2;0.03;0.02];%the main difference between two banks
Wassets2=[0.2;0.1;0.1;0.6];
%%
%Basel III requirements
%capital adequacy ratio >= 10.5% (minimum requirement+CCB) 
%RWA calculation cash, bond,loan to bank,loan to customer
Riskweights=[0;0;0.2;1;1];
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
% %%
% %optimization
% %CEO
% UCEO= CEO(Wassets,Rassets,Wliab,V,Eloss);%minus of utility function of CEO
% 
% fun = @(Wassets)CEO(Wassets,Rassets,Wliab,V,Eloss);
% 
% options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
%     'TolFun',1e-10,'MaxIter',1000000);  
% A=[A1;A2];
% B=[B1;B2];
% Aeq=ones(1,Nassets);
% Beq=1;
% LB=zeros(Nassets,1);
% UB=ones(Nassets,1);
% [optWCEO RESNORM MaxCEO exitflag] = fmincon(fun,Wassets,A,B,Aeq,Beq,LB,UB);
% MaxUCEO= -CEO(optWCEO,Rassets,Wliab,V,Eloss);
