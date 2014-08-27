clear all;
%Crank-Nicolson Method American Put
% Parameter Initialization
clear all;
Smax = 20; K =10; r = 0.04; sigma = 0.20; T =0.5;
deltat=0.002;
M=T/deltat;
deltax=sigma * sqrt(4*deltat);
N=ceil(log(Smax)/deltax);
Xgrid=(0:deltax:N*deltax)';
S=exp(Xgrid);
Nu=r-(sigma^2)*0.5;
%Crank-Nicolson
Pu = -.25 * deltat * (sigma^2/deltax^2 + Nu / deltax);
Pm = 1 + deltat * sigma^2/(2*deltax^2) + (r * deltat)*0.5;
Pd = -.25* deltat * (sigma^2/deltax^2 - Nu / deltax);
%Matrix generation
% Grid Generation
trans=eye(N+1)*Pm+diag(ones(N,1)*Pd,1)+diag(ones(N,1)*Pu,-1);
trans(1,1)=1;
trans(1,2)=-1;
trans(end,end-1)=1;
trans(end,end)=-1;
val=nan(N+1,M);
val(:,end)=max(0,K-S);
%Put
bottom = Xgrid(1) - Xgrid(2);
for m=M:(-1):2
  %Crank-Nicol-call/put
  z=-Pu*val(1:(end-2),m)-(Pm-2)*val(2:(end-1),m)-Pd*val(3:end,m);
  rhs=[0;z;bottom];
  val(:,m-1)=trans\rhs;
 %American
 val(:,m-1) = max(val(:,m-1), K-S);
 end
S0 = 10;
val(1:10, 1)
[C, P] = blsprice(S, K, r, T, sigma);
P(1:10)
