function utility= regulator(policy,X,Wequity,Rassets,X2,Wequity2,ASFratio,RSFratio,Covasset,Eloss,rho,Riskweights)
Wassets=X(1:5);
Wassets2=X2(1:5);
% Wliab=X(6:8,1);
% Wliab2=X2(6:8,1);
% Wequity=1-sum(Wliab,1);
% Wequity2=1-sum(Wliab2,1);


%Three constraints
RWAratio=policy(1);
NSFR=policy(2);
leverage=policy(3);%>=3%



%Optimization based on CEO's utility, bank 1
A1=Riskweights';
B1=Wequity/RWAratio;
A2=[1,1,1,1,1,0,0,0];
B2=Wequity/leverage;
% ASF=(Wequity+ASFratio*Wliab);
% RSF=(RSFratio*Wassets);%NSFR=ASF/RSF>=100%
A3=[RSFratio.*NSFR,-1.*ASFratio];
B3=Wequity;
A=[A1;A2;A3];
B=[B1;B2;B3];
Aeq=[1,1,1,1,1,0,0,0];;
Beq=Wequity;
LB=zeros(8,1);
UB=10000.*ones(8,1);
X
Wequity
UCEO= CEO(X,Wequity,Rassets,Covasset,Eloss);%minus of utility function of CEO

fun = @(X)CEO(X,Wequity,Rassets,Covasset,Eloss);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  
[optX RESNORM MaxCEO exitflag] = fmincon(fun,X,A,B,Aeq,Beq,LB,UB);

%Optimization based on CEO's utility, bank 2

A4=Riskweights';
B4=Wequity2/RWAratio;
A5=[1,1,1,1,1,0,0,0];
B5=Wequity2/leverage;
% ASF2=(Wequity2+ASFratio*Wliab2);
% RSF2=(RSFratio*Wassets2);
A6=[RSFratio.*NSFR,-1.*ASFratio];
B6=Wequity2;

A=[A4;A5;A6];
B=[B4;B5;B6];

Aeq2=[1,1,1,1,1,-1,-1,-1];
Beq2=Wequity;
LB2=zeros(8,1);
UB2=10000.*ones(8,1);

UCEO2= CEO(X2,Wequity2,Rassets,Covasset,Eloss);%minus of utility function of CEO

fun = @(X2)CEO(X2,Wequity2,Rassets,Covasset,Eloss);

options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',20000,...
    'TolFun',1e-10,'MaxIter',1000000);  
[optX2 RESNORM MaxCEO exitflag] = fmincon(fun,X2,A,B,Aeq,Beq,LB,UB);


optW=optX(1:5);
optW2=optX2(1:5);
%minimize rho*sigmaA*sigmaB

utility=(optW'*Covasset*optW)*(optW2'*Covasset*optW2);
optW
optW2
end