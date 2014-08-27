T=2;% Final Time as asked in the problem
N=5000;%No. of obs
R0=0.05;%Initial Value of X as given in the problem
delta=T/N;%creating delta for discretization of SDE
path=5000; % # of paths generated
sigma=0.08;
values=nan(path,1);%to contain values of each path
for j=1:path;%for loop ver no of path
z=randn(N,1);%generating random normal variables
R=nan(N,1);
R(1,:)=R0;
for i=1:N %for loop over time
    R(i+1,:)=0.85*(0.05-0.5*R(i,:))*delta+sigma*sqrt(delta)*z(i,:);
    %SDE as given in the problem
end
    EI=delta*sum(R);
end;
E1=mean(values);%finally computing the expected value of the function as asked
    %%%find 1/N*sum (exp 