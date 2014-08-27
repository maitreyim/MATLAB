function Project8
%Project8 Submission for Project8 in Computational Finance
%   Project 8 Questions and Solution
%   SUBMITTED BY ANSHUL MAHESHWARI

% %% ***********Q 1 a b c d e *************************
% clear all; clc;
% r0=0.05;si=0.1;kap=0.82;rmn=0.05;
% FV=1000;T=0.5;t=0;
% N=1000;% outer loop for pricing bond
% M=360;% time steps
% 
% %% **** Q 1 a ***
% BondP1a=VasBondP(r0,si,kap,rmn,FV,T,t,N,M);
% disp(['Q 1 a : Vasicek Value of Bond : ' num2str(BondP1a)]);
% 
% % %% **** Q 1 b ***
% clear T; clear M;
% C=30;T=[0.5:0.5:4]';
% [sizeT,~]=size(T);
% for i=1:sizeT-1
%     M = T(i)*366;
%     BondP(i)=VasBondP(r0,si,kap,rmn,C,T(i),t,N,M);%value of coupons at different times
% end
% M = T(sizeT)*366;
% BondP(sizeT)=VasBondP(r0,si,kap,rmn,C+FV,T(sizeT),t,N,M); %value of coupon + final payment
% CBondP1b=sum(BondP);
% disp(['Q 1 b : Vasicek Value of Coupon Paying Bond : ' num2str(CBondP1b)]);
% % % 
% % % %% ***** Q 1 c ***
% clear T;
% T=3/12;S=0.5;K=950;
% [EurCallOpt1c] = VasBondOptP(S,T,t,N,M,kap,rmn,si,r0,FV,K);
% disp(['Q 1 c : Vasicek Value of European Call Option pure discount Bond : ' num2str(EurCallOpt1c)]);
% 
% % % %% **** Q 1 d e 1st Method***
% clear T; clear M;clear S;
% C=30;N=100;
% Cou=C*ones(7,1);Cou(8,1)=1030;S=[0.5:0.5:4]';T=3/12;K=950;
% [CallpriceMet11d, CallpriceMet11e] = Q1deMet1VasBondOptP(S,T,t,N,kap,rmn,si,r0,FV,K,Cou);
% 
% disp(['Q 1 d : 1st Method Vasicek Value of European Call Option coupon paying Bond(Simulated) : ' num2str(CallpriceMet11d)]);
% disp(['Q 1 e : 1st Method Vasicek Value of European Call Option coupon paying Bond(Explicit) : ' num2str(CallpriceMet11e)]);
% 
% % % 
% % % %% **** Q 1 d e 2nd Method***
% clear T; clear M;clear S;
% C=30;
% Cou=C*ones(7,1);Cou(8,1)=1030;S=[0.5:0.5:4]';T=3/12;K=950;
% [Callprice1d, Callprice1e] = Q1deVasCPBOptPrice(S,kap,rmn,si,T,FV,Cou,C,N,K,t);
% 
% disp(['Q 1 d : 2nd Method Vasicek Value of European Call Option coupon paying Bond(Simulated) : ' num2str(Callprice1d)]);
% disp(['Q 1 e : 2nd Method Vasicek Value of European Call Option coupon paying Bond(Explicit) : ' num2str(Callprice1e)]);
% % %% ***** Q 2 a b c **********************************
% % % CIR MODEL Implementation
% % 
% % %%******* 2 a
% % % Calculate Option value using Euler's Integral
% % %*******
% clear all
% clear all; %clc;
% 
% % Initial Values
% r0=0.05;si=0.12;kap=0.92;rmn=0.055;
% t=0;K=950;T=0.5;FV=1000;S=1;
% NO=1000;% outer loop for pricing option
% MO=180;% time steps for options
% N=1000; % outer loop for Bond
% M=180; % time steps for Bond
% 
% EurCallOpt2a=CIRBondOptP(S,T,t,NO,MO,N,M,kap,rmn,si,r0,FV,K);
% disp(['Q 2 a : European call Option Value of Bond : ' num2str(EurCallOpt2a)]);
% 
% % % %%******** 2 b
% % % %% Calculate call Option using Implicit Finite Difference
% % % % %********
% % define the log price and time matrix
% Dt=0.002;Ds=0.001;DX=si*sqrt(3*Dt);CorP='C';
% [ Opt,Options ] = ImplicitFiniteDiff( r0,K,rmn,T,si,(r0-0.02), (r0+0.02),...
%     Ds,Dt,DX,CorP )
% 
% % % % 
% % % 
% % % %%******** 2 c 
% % % %% Compute price using Explicit formula
% % % 
% [CIROptPri2c] = CIRBondCallOptExplicit (t,T,kap,rmn,si,r0,S,K/1000,FV/1000 )*FV;
% disp(['Q 2 c : European call Option Value of Bond using Explicit : ' num2str(CIROptPri2c)]);

% % 
% %%******* Q 3 ********************************
% %% G2 ++ model implementation
clear all;%clc;

% Initial Values
x0=0;y0=0;phi0=0.03;r0=0.03;rho=0.7;a=0.1;b=0.3;si=0.03;eeta=0.08;phit=0.03;
t=0;K=1050;T=0.5;FV=1000;S=1;

NO=1000;% outer loop for pricing option
MO=366/2;% time steps for options
N=1000; % outer loop for Bond
M=366/2; % time steps for Bond

EurPutOpt3a=G2PPBondOptP(S,T,t,NO,MO,N,M,a,b,rho,eeta,phit,si,r0,FV,K,...
    x0,y0);
disp(['Q 3 : European Put Option Value of Bond (Simulated) : ' num2str(EurPutOpt3a)]);

% ************ Calc G2++ Explicit Price *****
[BondP3] = G2PPBondPExplicit(S,t,a,b,rho,eeta,phit,...
    si,x0,y0)*FV;
disp(['Q 3: Explicit Bond Price (Just Checking): ' num2str(BondP3)]);
[EuroPutOpt3b] = G2PPBondOptPExplicit(S,T,t,a,b,rho,eeta,phit,...
    si,x0,y0,K/1000)*FV;
disp(['Q 3: Explicit Bond Put Option Price (EXPLICIT) : ' num2str(EuroPutOpt3b)]);

%% END
end
% 
function [BondP] = CIRBondExplicit (t,T,kap,rmn,si,rt,FV )
% CIRBondExplicit CIR Model based Explicit price of Bond
% 

% calc h1, h2, h3
h1=sqrt(kap^2 + 2*si^2);h2=(kap+h1)/2;h3=2*kap*rmn/(si^2);

% calc B(t,T)
interExp=exp(h1*(T-t)) - 1;
B=interExp/(h2*interExp + h1);

%calc A(t,T)
num=h1*exp(h2*(T-t));
den=h2*interExp+h1;
A=(num/den)^h3;

%calculate bond price
BondP=FV*A*exp(-B*rt);
end

function [CIROptPri] = CIRBondCallOptExplicit (t,T,kap,rmn,si,rt,S,K,FV )
% CIRBondCallOptExplicit CIR Model based Explicit price of Call Option on Bond
% 

% calc h1, h2, h3
h1=sqrt(kap^2 + 2*si^2);h2=(kap+h1)/2;h3=2*kap*rmn/(si^2);

% calc B(T,S)
interExp=exp(h1*(S-T)) - 1;
BTS=interExp/(h2*interExp + h1);

%calc A(T,S)
num=h1*exp(h2*(S-T));
den=h2*interExp+h1;
ATS=(num/den)^h3;

%calc rstar
rstar=(log((ATS*FV)/K))/(BTS);

theta=sqrt(kap^2 + 2*si^2); 
phi=2*theta/(si^2*(exp(theta*(T-t))-1));
chi=(kap+theta)/si^2;

%calc PtS and PtT
PtS = CIRBondExplicit (t,S,kap,rmn,si,rt,FV );
PtT = CIRBondExplicit (t,T,kap,rmn,si,rt,FV );

% calc first chi-sq
P_chi_1=2*rstar*(phi+chi+BTS);
p1=4*kap*rmn/(si^2);
num=(2*phi^2)*rt*exp(theta*(T-t));
den=phi+chi+BTS;
q1=num/den;

%chi_1=ncx2cdf(P_chi_1,p1,q1);
chi_1=ncx2cdf(P_chi_1,p1,q1);
% 
% % calc second chi-sq
P_chi_2=2*rstar*(phi+chi);
p2=p1;
q2=num/(phi+chi);

chi_2=ncx2cdf(P_chi_2,p2,q2);

CIROptPri = PtS*chi_1 -K*PtT*chi_2;
end
% 
% 
function [BondP] = VasBondP(r0,si,kap,rmn,FV,T,t,N,M)
dt=(T-t)/M;
for i=1:N
    Z=randn(M,1);
    Vr(1)=r0;
    for j=2:M
        Vr(j)=Vr(j-1)+kap*(rmn-Vr(j-1))*dt+si*sqrt(dt)*Z(j-1);
    end
    VI(i)=dt*sum(Vr);
end
BondP=(1/N)*sum(exp(-VI))*FV;
end
% 
function [BondP] = VasBondPExplicit(r0,si,kap,rmn,FV,T,t)

BtT=(1/kap)*(1-exp(-kap*(T-t)));
AtT=exp((rmn-0.5*(si/kap)^2)*(BtT - (T-t)) - ((si^2)/4*kap)*(BtT^2));

BondP=AtT*exp(-BtT*r0)*FV;
end

function [EurCallOpt] = VasBondOptP(S,T,t,NO,MO,kap,rmn,si,r0,FV,K)
% Calculate path of rate till time T
dtO=(T-t)/MO;
for i=1:NO
    ZO=randn(MO,1);
    VrO(1)=r0;
    for j=2:MO
        VrO(j)=VrO(j-1)+kap*(rmn-VrO(j-1))*dtO+si*sqrt(VrO(j-1))*sqrt(dtO)*ZO(j-1);
    end
    % for each end rs calc N number of Bond Price
    BondP(i)=VasBondP(VrO(MO),si,kap,rmn,FV,S,T,NO,MO);
    VIO(i)=dtO*sum(VrO);
end
EurCallOpt=(1/NO)*sum((exp(-VIO)).*(max(BondP-K,0)));
end
%%
function [EurCallOpt] = CIRBondOptP(S,T,t,NO,MO,N,M,kap,rmn,si,r0,FV,K)
% % Calculate path of rate till time T
dtO=(T-t)/MO; dt=(S-T)/M;
for i=1:NO
    ZO=randn(MO,1);
    VrO(1)=r0;
    for j=2:MO
        VrO(j)=VrO(j-1)+kap*(rmn-VrO(j-1))*dtO+si*sqrt(VrO(j-1))*sqrt(dtO)*ZO(j-1);
    end
    % for each end rs calc N number of Bond Price
    for j=1:N
        Z=randn(M,1);
        Vr(1)=VrO(MO);
        for k=2:M
            Vr(k)=Vr(k-1)+kap*(rmn-Vr(k-1))*dt+si*sqrt(Vr(k-1))*sqrt(dt)*Z(k-1);
        end
        VI(j)=dt*sum(Vr);
    end
    BondP(i)=(1/N)*sum(exp(-VI))*FV;
    VIO(i)=dtO*sum(VrO);
end
EurCallOpt=(1/NO)*sum((exp(-VIO)).*(max(BondP-K,0)));
end
% %%
function [EurPutOpt] = G2PPBondOptP(S,T,t,NO,MO,N,M,a,b,rho,eeta,phit,...
    si,r0,FV,K,x0,y0)
% Calculate path of rate till time T
dtO=(T-t)/MO; dt=(S-T)/M;
for i=1:NO
    ZOX=randn(MO,1); ZO2=randn(MO,1);
    ZOY=rho.*ZOX + sqrt(1-rho^2).*ZO2;
    VxO(1)=x0;VyO(1)=y0;VrO(1)=r0;
    for j=2:MO
        VxO(j)=VxO(j-1)+ (-1)*a*VxO(j-1)*dtO + si*sqrt(dtO)*ZOX(j-1);
        VyO(j)=VyO(j-1)+ (-1)*b*VyO(j-1)*dtO + eeta*sqrt(dtO)*ZOY(j-1);
        VrO(j)=phit+ VxO(j)+VyO(j);
    end
    % for each end rs calc N number of Bond Price
    for j=1:N
        ZX=randn(M,1); Z2=randn(M,1);
        ZY=rho.*ZX + sqrt(1-rho^2).*Z2;
        Vx(1)=VxO(MO);Vy(1)=VyO(MO);Vr(1)=VrO(MO);%+phit;
        for k=2:M
            Vx(k)=Vx(k-1)+ (-1)*a*Vx(k-1)*dt + si*sqrt(dt)*ZX(k-1);
            Vy(k)=Vy(k-1)+ (-1)*b*Vy(k-1)*dt + eeta*sqrt(dt)*ZY(k-1);
            Vr(k)=phit+Vx(k)+Vy(k);
        end
        VI(j)=dt*sum(Vr);
    end
    BondP(i)=(1/N)*sum(exp(-VI))*FV;
    VIO(i)=dtO*sum(VrO);
end
EurPutOpt=(1/NO)*sum((exp(-VIO)).*(max(K-BondP,0)));
%EurPutOpt=(1/NO)*sum((exp(-VIO)).*(max(BondP-K,0)));
end
% %%
function [BondP] = G2PPBondPExplicit(T,t,a,b,rho,eeta,phit,...
    si,x0,y0)
% 
% % calc V(t,T)
first = ((si/a)^2)*(T - t + (2/a)*exp(-a*(T-t)) - (1/(2*a))*exp(-2*a*(T-t))...
    - 3/(2*a));
second = ((eeta/b)^2)*(T - t + (2/b)*exp(-b*(T-t)) - (1/(2*b))*exp(-2*b*(T-t))...
    - 3/(2*b));
third= (2*rho*si*eeta/(a*b)) * (T - t + (exp(-a*(T-t))-1)/a ...
    + (exp(-b*(T-t)) - 1)/b - (exp(-(a+b)*(T-t))-1)/(a+b));

VtT = first + second + third;
% 
% % Calculate P(t,T)
% % since phit is constant thorughout the integral will be phit*(T-t)
fi = phit*(T-t);
se = (1-exp(-a*(T-t)))*x0/a;
thi = (1-exp(-b*(T-t)))*y0/b;
fo = 0.5*VtT;

BondP = exp(-fi - se - thi + fo);
end
%%
function [EuroPutOpt] = G2PPBondOptPExplicit(S,T,t,a,b,rho,eeta,phit,...
    si,x0,y0,K)
% 
% % find Big sigma
first = (0.5*(si^2)/(a^3)) * ((1 - exp(-a*(S-T)))^2) * (1 - exp(-2*a*(T-t)));
second = (0.5*(eeta^2)/(b^3)) * ((1 - exp(-b*(S-T)))^2) * ...
    (1 - exp(-2*b*(T-t)));
third = (2*rho*si*eeta/(a*b*(a+b)))*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))...
    *(1-exp(-(a+b)*(T-t)));

sigma=sqrt(first + second + third);
% 
% % find Put Option Price
% get P(t,S)
PtS = G2PPBondPExplicit(S,t,a,b,rho,eeta,phit,...
    si,x0,y0);
% get P(t,T)
PtT = G2PPBondPExplicit(T,t,a,b,rho,eeta,phit,...
    si,x0,y0);

intTerm = log(K*PtT/PtS)/sigma;
Nd1 = N(intTerm - 0.5*sigma); Nd2=N(intTerm + 0.5*sigma);
Nd1 = normcdf(intTerm-0.5*sigma,0,1); Nd2=normcdf(intTerm+0.5*sigma,0,1);

EuroPutOpt=-PtS*Nd1 + PtT*K*Nd2;
end
% %%
function [Callprice1d, Callprice1e] = Q1deVasCPBOptPrice(S,kap,rmn,si,T,FV,Cou,C,N,K,t)

dt=S-T;
%dt=flipud(dt);% flip the matrix upside down
[sizeS,~]=size(S);
% Calculate As and Bs for each coupon payment day

B= (1/kap).*(1-exp(-kap.*dt));
B2=B.^2;
A= exp((rmn- (si^2)/(2*kap^2)).*(B- dt)- ((si^2)/(4*kap)).*B2);

% solve for r* (R)
syms R
[R]=solve(C*(A(1)*exp(-B(1)*R)+...
A(2)*exp(-B(2)*R)+...
A(3)*exp(-B(3)*R)+...
A(4)*exp(-B(4)*R)+...
A(5)*exp(-B(5)*R)+...
A(6)*exp(-B(6)*R)+...
A(7)*exp(-B(7)*R)+...
A(8)*exp(-B(8)*R))+...
FV*(A(8)*exp(-B(8)*R))...
-K);

% get Strike price of $ 1 for each coupon payment day and also find price
% of $1 bonds at different maturities
KC= zeros(sizeS,1); % strike prices - preallocating matrix
BondP1= zeros(sizeS,1); % Bond prices - preallocating matrix
BondP1Exp= zeros(sizeS,1); % Bond prices - preallocating matrix
R = double(R);
KC= double(A.*exp(-B.*R));
M = S.*360;
for i=1:sizeS
    BondP1(i)=VasBondP(R,si,kap,rmn,1,S(i),t,N,M(i));%value of coupons+Principal at different times using Simulation
    BondP1Exp(i) = VasBondPExplicit(R,si,kap,rmn,1,S(i),t);%value of coupons+Principal at different times using Explicit
end
MT=T.*360;
BondPT=VasBondP(R,si,kap,rmn,1,T,t,N,MT);
BondPTExp = VasBondPExplicit(R,si,kap,rmn,1,T,t);
% get strike price based on coupon amount   
KOpt= KC.*(Cou); % multiply coupon payments to get Final strike price
% 
% calculate sigmap, d1, d2 
sip= si.*((1-exp(-kap.*(dt)))/kap).*sqrt((1-exp(-2.*kap.*dt))/(2*kap));
d1= log(BondP1./((KOpt).*BondPT))./(sip) +sip/2;%simulated bond
d2= d1-sip;%simulated Bond

dExp1= log(BondP1Exp./((KOpt).*BondPTExp))./(sip) +sip/2;%explicit bond
dExp2= d1-sip;%explicit Bond

% calculate option prices
c=zeros(sizeS,1);cExp=zeros(sizeS,1);
for i=1:sizeS
    c(i)=max(BondP1(i).*normcdf(d1(i),0,1)- KOpt(i).*...
        BondPT.*normcdf(d2(i),0,1),0); %vector of option price at each coupon day for simulated bond
    cExp(i)=max(BondP1Exp(i).*normcdf(dExp1(i),0,1)- KOpt(i).*...
        BondPTExp.*normcdf(dExp2(i),0,1),0); %vector of option price at each coupon day for explicit bond
end
Callprice1d= sum(c.*Cou);% multiply coupons with option price. SUM to arrive at call price.Simulated
Callprice1e= sum(cExp.*Cou);% multiply coupons with option price. SUM to arrive at call price. Explicit
end
% 
% %%
function [CallpriceMet11d, CallpriceMet11e] = Q1deMet1VasBondOptP(S,T,t,N,kap,rmn,si,r0,FV,K,Cou)
% Calculate path of rate till time T
MO=(T-t)*360;
dtO=(T-t)/MO;
[sizeS,~]=size(S);
M=(S-T)*360;
for i=1:N
    ZO=randn(MO,1);
    VrO(1)=r0;
    for j=2:MO
        VrO(j)=VrO(j-1)+kap*(rmn-VrO(j-1))*dtO+si*sqrt(VrO(j-1))*sqrt(dtO)*ZO(j-1);
    end
    % for each end rs calc N number of Bond Price
    for k=1:sizeS
        BondP(k)=VasBondP(VrO(MO),si,kap,rmn,Cou(k),S(k),T,N,M(k));
        BondPExp(k)=VasBondPExplicit(VrO(MO),si,kap,rmn,Cou(k),S(k),T);%value of coupons+Principal at different times using Explicit
    end
    BondFullP(i)=sum(BondP);
    BondFullPExp(i)=sum(BondPExp);
    VIO(i)=dtO*sum(VrO);
end
CallpriceMet11d=(1/N)*sum((exp(-VIO)).*(max(BondFullP-K,0)));
CallpriceMet11e=(1/N)*sum((exp(-VIO)).*(max(BondFullPExp-K,0)));
end
% 
% 
% 
% 
% 
% 
% 






