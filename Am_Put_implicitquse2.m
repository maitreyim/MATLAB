clear all;
% Implicit Finite Difference Method for American Put Option
% Parameter Initialization
Smax = 75; K =10; r = 0.04; sigma = 0.20; T =0.5;
deltat=0.002;
M=T/deltat;
% deltax=.1;
deltax=sigma * sqrt(3*deltat);
% N=ceil(Smax/deltax);
N=500;
Sgrid=linspace(Smax,0,N)';
Nu=r-(sigma^2)*0.5;
%Implicit
% Pu= deltat*(sigma^2/(2*deltax^2)+Nu/(2*deltax));
% Pm= 1-deltat*sigma^2/(deltax^2)-r*deltat;
% Pd= deltat*(sigma^2/(2*deltax^2)-Nu/(2*deltax));
%Explicit
% Pu = -.5 * deltat * (sigma^2/deltax^2 + Nu / deltax);
% Pm = 1 + deltat * sigma^2/deltax^2 + r * deltat;
% Pd = -.5 * deltat * (sigma^2/deltax^2 - Nu / deltax);
%Crank-Nicolson
Pu = -.25 * deltat * (sigma^2/deltax^2 + Nu / deltax);
Pm = 1 + deltat * sigma^2/(2*deltax^2) + (r * deltat)*0.5;
Pd = -.25* deltat * (sigma^2/deltax^2 - Nu / deltax);
trans=eye(N)*Pm+diag(ones(N-1,1)*Pd,1)+diag(ones(N-1,1)*Pu,-1);
trans(1,1)=1;
trans(1,2)=-1;
trans(end,end-1)=1;
trans(end,end)=-1;
val=nan(N,M);
% val(:,end)=max(0,Sgrid-K);
val(:,end)=max(0,K-Sgrid);
top=Sgrid(1)-Sgrid(2);
%Put
bottom = Sgrid(end) - Sgrid(end-1);
%For European Put
for m=M:(-1):2
    %Call-imp/ex
%    rhs=[top;val(2:(end-1),m);0];
%Put-Implicit/Explicit
%     rhs=[0;val(2:(end-1),m);bottom];
%Crank-Nicol-call/put
    z=-Pu*val(1:(end-2),m)-(Pm-2)*val(2:(end-1),m)-Pd*val(3:end,m);
    %put-eu/am
    rhs=[0;z;bottom];
    val(:,m-1)=trans\rhs;
    %American
    %val(:,m-1) = max(val(:,m-1), K-Sgrid);
      
end
S0 = 11;
[toss, index] = min(abs(S0 - Sgrid));
val(index, 1)
[C, P] = blsprice(S0, K, r, T, sigma);










