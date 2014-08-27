%Question6: Calculate The Integral
%Part A: Euler's Discretization Method
n=1000;
alpha=0.5;
for i=1:n
    f(i)=sqrt(1-(i/n)^2);
end
z=(1/n)*sum(f);
I=4*z;
%Part B: Calculate the Integral Using Monte Carlo
% A:Generate Uniformly Distributed Random Numbers using LGM Algorithm
n = 1000;%number asked for in the question
x = nan(n, 1);% pre allocating memory space for x's
u = nan(n, 1);% pre allocating memory space for u's
x(1) = 10; % seed/staring number for x's
m = 2^31 - 1;%standard m for LGM algorithm
b=0;%standard b for LGM algorithm
a= 7^5; % standard a in LGM algorithm
for i=1:n
x(i+1) = mod((a*x(i)),m);
%Generate U(i)
u(i)= x(i)/m;
s(i)=sqrt(1-u(i)^2);
%Part B: Calculate the Integral using Importance Sampling Method
h(i)=(1-alpha*(u(i)^2))/(1-alpha/3);
g(i)=sqrt(1-u(i)^2);
k(i)=g(i)/h(i);
end
s1=(1/n)*sum(s);
I1=4*s1;
k1=(1/n)*sum(k);
