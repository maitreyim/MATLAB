% Problem 2, Part A
%Compute E(1+X3)^(1/3)
T=3; % Final Time value as given
N=1000;% # of obs
rng(1);
z=sqrt(3)*randn(N,1);% First set of standard wiener process
rng(2);
w=sqrt(3)*randn(N,1);% Second set of standard wiener process
a=1/4;% value of a as given in the problem
b=1/3;% value of b as given in the problem
path=5000; % # of paths
delta=T/N;
X0=1; % initial value as given in the problem
for j=1:path;   % loop over path
z=randn(N,1);
X = zeros(N,1);
X(1)=X0;
for i=1:N % loop over time
    X(i+1,:)=X(i,:)+a*X(i,:)*(delta)+b*X(i,:)*sqrt(delta)*z(i,:)-(3/4)*X(i,:)*sqrt(delta)*z(i,:);
    %SDE as given
        end
    values(j,1)=X(N,1);
end;
E3=mean((1+values).^(1/3));% computing the expected value

