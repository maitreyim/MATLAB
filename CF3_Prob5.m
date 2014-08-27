%Problem 5
%Part A
% A:Generate Uniformly Distributed Random Numbers using LGM Algorithm
n = 100;%number asked for in the question
x = nan(2*n, 1);% pre allocating memory space for x's
u = nan(2*n, 1);% pre allocating memory space for u's
x(1) = 10; % seed/staring number for x's
m = 2^31 - 1;%standard m for LGM algorithm
b=0;%standard b for LGM algorithm
a= 7^5; % standard a in LGM algorithm
for i=1:(2*n)
x(i+1) = mod((a*x(i)),m);
%Generate U(i)
u(i)= x(i)/m;
end
z=reshape(u, n, 2);
scatter(z,n,2);
%plot(z);


