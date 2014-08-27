clear all; clc;
r0=0.05;si=0.1;kap=0.82;rmn=0.05;
FV=1000;T=0.5;t=0;
N=1000;% outer loop for pricing bond
M=360;% time steps

%% **** Q 1 a ***
BondP1a=VasBondP(r0,si,kap,rmn,FV,T,t,N,M);
disp(['Q 1 a : Vasicek Value of Bond : ' num2str(BondP1a)]);
% 
% %% **** Q 1 b ***
clear T; clear M;
C=30;T=[0.5:0.5:4]';
[sizeT,~]=size(T);
for i=1:sizeT-1
    M = T(i)*366;
    BondP(i)=VasBondP(r0,si,kap,rmn,C,T(i),t,N,M);%value of coupons at different times
end
M = T(sizeT)*366;
BondP(sizeT)=VasBondP(r0,si,kap,rmn,C+FV,T(sizeT),t,N,M); %value of coupon + final payment
CBondP1b=sum(BondP);
disp(['Q 1 b : Vasicek Value of Coupon Paying Bond : ' num2str(CBondP1b)]);
% 
% %% ***** Q 1 c ***
clear T;
T=3/12;S=0.5;K=950;
[EurCallOpt1c] = VasBondOptP(S,T,t,N,M,kap,rmn,si,r0,FV,K);
disp(['Q 1 c : Vasicek Value of European Call Option pure discount Bond : ' num2str(EurCallOpt1c)]);