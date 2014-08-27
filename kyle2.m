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
  X(i+1,1) = X(i,1) + (1/5 - X(i,1)/2) * delta + (2/3)*sqrt(delta)*z(i,1);
end
    values(j,1)=X(N,1)^(1/3 );
end;
E1=mean(values);
    