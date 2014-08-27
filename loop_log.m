nSteps = [10, 20, 40, 80, 100, 200, 500]';% Different Time Steps as specified in the problem
vals = nan(size(nSteps)); % Pre allocating memory space for Different European Option Values obtained through Trinomial Pricing
S=32;
r=0.05;
sigma=0.24;
K=30;
T=0.5;
for i= 1:size(nSteps,1)
    
  %Get the price of a call option using Trinomial tree and log prices
   vals(i) = CF4_Trino_log(log(S), log(K), r, T, sigma, nSteps(i));
end;
%Plot trinomial 
figure(1)
plot(nSteps,vals);
Title('Trinomial w.r.t # of steps');
Xlabel('# of steps');


