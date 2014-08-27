T=2;
N=1000;
X0=1;
delta=T/N;
path=50;
values=nan(path,1);
for j=1:path;
z=randn(N,1);
X=nan(N,1);
X(1,:)=X0;
for i=1:N
    X(i+1,:)=X(i,:)+((1/3)*X(i,:)^2*(1/5 - 1/2*X(i,:)^3)+(.5*(4/9)*-2/9)*X(i,:)^(-5))*delta + ...
       ((2/9)*X(i,:)^(-2) )*sqrt(delta)*z(i,:);
end
    values(j,1)=X(N,1);
end;
E1=mean(values);
    