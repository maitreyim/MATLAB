S0_set = 1:1:20;,
K=10,
r=0.04;
T=0.5;
sigma=0.2;
for i = 1:size(S0_set,2)
    
vals(i) = CF4_AP(S0_set(i),T,r,sigma,K );
end
vals