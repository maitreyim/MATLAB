clear all;
% Implicit Finite Difference Method for American Put Option
% Parameter Initialization
Smax = 20; K =10; r = 0.04; sigma = 0.20; T =0.5;
deltat=0.002;
M=T/deltat;
%deltax=1;
deltax=sigma * sqrt(4*deltat);
N=ceil(log(Smax)/deltax);
%Sgrid=flipud(linspace(log(Smax),0,N)');
Xgrid=(0:deltax:N*deltax)';
S=exp(Xgrid);
Nu=r-(sigma^2)*0.5;
%Implicit
Pu= -deltat*(sigma^2/(2*deltax^2)+Nu/(2*deltax));
Pm= 1+deltat*sigma^2/(deltax^2)+r*deltat;
Pd= -deltat*(sigma^2/(2*deltax^2)-Nu/(2*deltax));
% Grid Generation
trans=eye(N+1)*Pm+diag(ones(N,1)*Pd,1)+diag(ones(N,1)*Pu,-1);
trans(1,1)=1;
trans(1,2)=-1;
trans(end,end-1)=1;
trans(end,end)=-1;
val=nan(N+1,M);
val(:,end)=max(0,K-S);
%Put
%bottom = Sgrid(end) - Sgrid(end-1);
bottom = Xgrid(1) - Xgrid(2);
for m=M:(-1):2
%Put-imp/ex
 rhs=[0;val((end-1):-1:2,m);bottom];
 val(:,m-1)=flipud(trans\rhs);
%American
 val(:,m-1) = max(val(:,m-1), K-S);
end
S0 = 10;
% [toss, index] = min(abs(S0 - S));
val(1:10, 1)
[C, P] = blsprice(S, K, r, T, sigma);
P(1:10)
