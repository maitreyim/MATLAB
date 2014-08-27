%Homework Assignment #1 For Computational Methods in Finance Class
%Problem #1
% A:Generate Uniformly Distributed Random Numbers using LGM Algorithm
n = 10000;%number asked for in the question
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
end
j=20;%defining bins
figure(1)
hist(u,j);% drawing the histogram of u's
xlabel('Bin','Fontsize',12)
ylabel('Frequency','Fontsize',12)
Title('Histogram of Uniformly Generated Random Numbers Using LGM ALgorithm','Fontsize',12)
% B:Generate Uniformly Distributed Random Numbers using RAND of MATLAB
r=rand(10000,1);% RAND function
j=20;%defining bins
figure(2)
hist(r,j);% drawing histogram
xlabel('Bin','Fontsize',12)
ylabel('Frequency','Fontsize',12)
Title('Histogram of Uniformly Generated Random Numbers Using RAND Function','Fontsize',12)
% Compare the two distributions as generated above by Kolmogorov-Smirnoff
% Test
%Two Sample KS Test
[h,p,k]=kstest2(u,r);
% Plot CDFs of two above generated distributions
figure(3)
CDF1 = cdfplot(u);
hold on
CDF2 = cdfplot(r);
set(CDF1,'Color','r');
set(CDF2,'Color','b');
legend([CDF1 CDF2],'CDF1(Uniform-LGM)','CDF2(Uniform-RAND)','Location','SE')
Title('Comparison of Two Uniformly Distributed Random Numbers','Fontsize',12)
