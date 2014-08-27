% Problem 1, Part c
% Compute Expected value of (X2Y2 1(X2>1))
T=2;%Final Time as asked in the question
N=1000;%# of obs
X0=1;%Initial Value of X
Y0=3/4;%Initial Value of Y
delta=T/N;% Creating Delta
path=500;%# of paths generated
values=nan(path,1);%pre allocating memory space
for j=1:path;% loop over path
z=randn(N,1);% creating first sequence of random normal numbers of Zt
w=randn(N,1);% Creating second sequence of random normal numbers of Wt
X=nan(N,1);
Y=nan(N,1);
X(1)=X0;
Y(1)=Y0;
for i=1:N % time loop
    X(i+1,:)=X(i,:)+(1/5-0.5*X(i,:))*delta+(2/3)*sqrt(delta)*z(i,:);% As asked in the question
     Y(i+1,:)=Y(i,:)+((2/(1+(i-1)*delta))*Y(i,1)+ (1+(i-1)*delta^3)/3)*delta+(1+((i-1)*(delta^3))/3)*sqrt(delta)*z(i,1);
     %as asked in the question
end
    
values(j,:)=X(N,1)*Y(N,1)*(X(N,1)>1);
end;
P=mean(values);% Finally computing the expected value

   