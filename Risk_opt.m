clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%Assumptions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%OPTIMIZED PORTFOLIO%%%%%
Wassets=[ 0.425257732;	0.180412371;	0.077319588;	0.103092784;	0.056701031;	0.06443299;	0.067010309;	0.025773196];
 % Importing Data 
[NUMERIC,TXT,RAW]=xlsread('C:\Users\maitreyi\Desktop\risk_opt.xlsx','Data','b2:i316');
Assetret=cell2mat(RAW);
return1=Assetret(:,1);
return2=Assetret(:,2);
return3=Assetret(:,3);
return4=Assetret(:,4);
return5=Assetret(:,5);
return6=Assetret(:,6);
return7=Assetret(:,7);
return8=Assetret(:,8);
%Covariance
c= cov(Assetret);
stdev=std(Assetret);
m=mean(Assetret);
[PortRisk, PortReturn, PortWts] = portopt(m, c);
scatter(stdev, m)
hold on
plot(PortRisk,PortReturn,'DisplayName','PortReturn vs. PortRisk','XDataSource','PortRisk','YDataSource','PortReturn');figure(gcf)
xlabel('Risk (Standard Deviation)')
ylabel('Expected Return')
title('Mean-Variance-Efficient Frontier')
grid on
