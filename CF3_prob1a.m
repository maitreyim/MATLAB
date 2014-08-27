%Question 1, Part A
%Evaluate the expected value of X2^(1/3)
T=2;% Final Time as asked in the problem
N=5000;%No. of obs
X0=1;%Initial Value of X as given in the problem
delta=T/N;%creating delta for discretization of SDE
path=5000; % # of paths generated
values=nan(path,1);%to contain values of each path
for j=1:path;%for loop ver no of path
z=randn(N,1);%generating random normal variables
X=nan(N,1);
X(1,:)=X0;
for i=1:N %for loop over time
    X(i+1,:)=X(i,:)+((1/5)-0.5*X(i,:))*delta+(2/3)*sqrt(delta)*z(i,:);
    %SDE as given in the problem
end
    values(j,1)=nthroot (X(N,1),-3);%taking only real root
end;
E1=mean(values);%finally computing the expected value of the function as asked
    