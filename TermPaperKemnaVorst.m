%PathAvgOption Submission for Term Paper Initial function to replicate
%   Table 1 results from Kemna-Vorst paper of 1989
% Average path options - Asian Options
%% Define Headers of Table
Price{1,1}='Table 1';
Price{2,1}='Monte Carlo simulation results';
ctr=3;
Price{ctr,1}={'sigma'}; Price{ctr,2}={'K'};
Price{ctr,3}={'Ca'}; Price{ctr,4}={'C(rnd)b'}; Price{ctr,6}={'C(rnd) red. var c'};
%% Define input parameters
r=[.03 .05 .07]; si=[.2 .3 .4]; K=[35 40 45]; T=4/12; S0=40;%Cox and Rubinstein 1985 paper p. 216 values
%% Define Monte Carlo parameters
N=10000;%number of simulations
n=T*12*22-1;%number of time steps based on trading days in T*12 months
%% build the table
for i=1:3 % rate loop
    ctr=ctr+1;
    Price{ctr,1}={'r='};Price{ctr,2}= 1+r(i);
    for j=1:3 % volatility/sigma loop
        for k=1:3 % Strike Price Loop
            ctr=ctr+1;
            Price{ctr,1}=si(j);
            Price{ctr,2}=K(k);
            [Ca Cb siCb CcRedVar siCc]=getPathAvgOptPrc(N,n,T,S0,r(i),si(j),K(k));
            Price{ctr,3}=num2str(Ca); 
            Price{ctr,4}=num2str(Cb); Price{ctr,5}= ['"(' num2str(siCb) ')'];
            Price{ctr,6}=num2str(CcRedVar); Price{ctr,7}=['"(' num2str(siCc) ')'];
        end
    end
end
Price %print Table 1 results