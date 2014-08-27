%=========================================================================
% Solution for BDT Model
% Priyank Gandhi
%=========================================================================
clear; clc; close all;
global idx1;
global bdt_tree;

% Load all the required files
prices = xlsread('C:\Users\maitreyi\Desktop\Fixed_Income\Hw6_pfilea',.....
    'sheet1','A1:A29');
vol = xlsread('C:\Users\maitreyi\Desktop\Fixed_Income\Hw6_voldat',....
'sheet1','A1:A29');

% Load an empty file for the tree. You will need a T x T matrix. In this
% case we have given you data out to 15 years in steps of 0.5 years.
% Therefore you should have a 30 x 30 matrix for the tree
bdt_tree = zeros(size(prices,1), size(prices,1));

% For the first value you do not need to do an optimization. This is easily
% calculated by setting D(0.5) = 1/(1+r/2)
bdt_tree(1,1) = 2*(1/prices(1) - 1);

% Set the parameters for the optimization
LowerBound = (10^(-9)) .* ones(5,1);
UpperBound = 1000000 .* ones(5,1);
Opti = optimset('fzero');
Opti.MaxIter = 1000;
Opti.TolFun = 10^(-5);
Opti.TolX = 10^(-8);
Opti.Display = 'off';
Opti.MaxFunEvals = 1000000;

%Set the initial value of the interest rate
r_0 = 0.05;

% Now run a loop for the remaining 29 maturities
for idx1 = 2 : size(prices,1)
        
    % Run the optimization.Use the fzero function so that you match the
    % D(T) function of each maturity exactly. See the objective function
    % objective_bdt
    
    r = fzero('objective_bdt', r_0, Opti);
    
    % Use the optimized value to fill the tree
    bdt_tree(1, idx1) = r;
    for idx2 = 2 : idx1
        bdt_tree(idx2, idx1) = bdt_tree(1, idx1)*exp(-2*(idx2-1)*vol(idx1-1)*sqrt(0.5));
    end    
    
    % Return to the next optimization
    
end

% At the end of this process you should have the bdt tree. Look at the bdt
% tree to make sure numbers are close to those given to you. And that the
% tree has not exploded with increased maturity. A tree with very high
% interest rates usually means that you have made an error.

% Now it is easy to compute the expected value and compare that to the 
% forward rates