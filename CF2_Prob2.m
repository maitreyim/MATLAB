%Problem 2: Calculate the Expected value of (x^2+sin(y)+x^2y) by simulation
rng(1);
N=5000000;
z1=randn(N,1);%Generating first series of normal distribution
rng(2)
z2=randn(N,1);%Generating second series of normal distribution
x=z1;%x=mean1+sigma1*z1,mean1=0,sigma1=1,so x=z1.
y=0.65*z1+sqrt(1-(0.65^2))*z2;%y=mean2+sigma2*rho*z1+sigma2*sqrt(1-rho^2)*z2
%mean2=0,sigma2=1,rho=0.65 from the given var-covar matrix.
a=mean((x.^3+(x.^2).*y+sin(y)));
