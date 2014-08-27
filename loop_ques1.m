clear all
nSteps = [10,15, 20, 40, 80, 100, 200, 500]';% Different Time Steps as specified in the problem
vals = nan(size(nSteps)); % Pre allocating memory space for Different European Option Values obtained through Trinomial Pricing
S0=32;
r=0.05;
sigma=0.24;
K=30;
T=0.5;
Acc_val=blsprice(32, 30,0.05, 0.6, 0.24);

for i = 1:size(nSteps,1)

vals(i)=bino_diff2(S0,sigma,nSteps(i),r,K );
%error(i)=Acc_val-vals(i);

end
plot(nSteps,vals);
%Legend('European Call Price dynamics with # of steps');
Xlabel('# of steps');
Ylabel('European Call Price Dynamics 2');
