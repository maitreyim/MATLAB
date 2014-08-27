1;

clear all;
%Demonstration of effectiveness of control variate technique
%in reducing variance of arithmetic asian option price

%%%%% Input parameters %%%%%%%%
S=100; %spot price
K=105; %strike
T=1; %maturity
r=0.05; %interest rate
sigma=0.25; %volatility
NSimulations=100000; %no of monte carlo simulations
NSteps=100; %no of time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=T/(NSteps-1);
vsqrdt=sigma*dt^0.5;
drift=(r-(sigma^2)/2)*dt;
x=randn(NSimulations,NSteps);
Smat=zeros(NSimulations,NSteps);
Smat(:,1)=S;
for i=2:NSteps,
	Smat(:,i)=Smat(:,i-1).*exp(drift+vsqrdt*x(:,i));
end

%TerminalSVec=Smat(:,NSteps);
AverageSVec=mean(Smat,2);
printf("[LowerLimit MCPrice UpperLimit]\n");

%*************  calculate call price ***********
callpayoffvec=max(AverageSVec-K,0);
mccallprice=exp(-r*T)*mean(callpayoffvec);
std_dev_call=(var(callpayoffvec))^0.5;
%calculate bounds of expected value assuming 95% confidence level
mclimit_call=1.96*std_dev_call/NSimulations^0.5;
call_lowerlimit=mccallprice-mclimit_call;
call_upperlimit=mccallprice+mclimit_call;

printf("Call Prices : [%f %f %f]\n",call_lowerlimit,mccallprice,call_upperlimit);

%*************  calculate put price ***********

putpayoffvec=max(K-AverageSVec,0);
mcputprice=exp(-r*T)*mean(putpayoffvec);
std_dev_put=(var(putpayoffvec))^0.5;
%calculate bounds of expected value assuming 95% confidence level
mclimit_put=1.96*std_dev_put/NSimulations^0.5;
put_lowerlimit=mcputprice-mclimit_put;
put_upperlimit=mcputprice+mclimit_put;

printf("Put Prices : [%f %f %f]\n",put_lowerlimit,mcputprice,put_upperlimit);







