%Problem1: Calculate the correlation coefficient by Monte Carlo Simulation
rng(1);
N=5000000;
z1=randn(N,1);%Generating First Series of Normal Distribution
rng(2);
z2=randn(N,1);%Generating Second Series of Normal Distribution
x=z1;%x=mean1+sigma1*z1,mean1=0,sigma1=1,so x=z1.
y=0.7*z1+sqrt(1-(0.7^2))*z2;%y=mean2+sigma2*rho*z1+sigma2*sqrt(1-rho^2)*z2
%mean2=0,sigma2=1,rho=0.7 from the given var-covar matrix.
xMean=mean(x);
yMean=mean(y);
for i=1:N
    numerator(i)=(x(i)-xMean)*(y(i)-yMean);
    xDenom(i) = (x(i)-xMean)*(x(i)-xMean);
    yDenom(i) = (y(i)-yMean)*(y(i)-yMean);
end;
rho=(sum(numerator)/N)/(sqrt(sum(xDenom)/N)*sqrt(sum(yDenom)/N));