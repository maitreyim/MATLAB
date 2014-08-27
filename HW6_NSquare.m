% HW 6 :: Submitted by Jie "Edward", Anshul, Sally, Chao
% Calculation of N^2 values of interest rates
% MAIN Script to run off the calculations

% clean up
clear all; clc;

% STEP 1:: Choose your alpha, beta and sigmas
alphax = .0052; betax = 0.053; sigmax= 0.0088;
% alphay = 0; don't need alphay as it is assumed to be 0
betay = 0.65; sigmay=0.024;

% read data from excel file
cmt = xlsread('Homework 5 data.xlsx','K2:P651')/100;

x0=[alphax, betax, sigmax, betay, sigmay];
f = @(x)calcVasicekNSquareRMSE(x,cmt);

my_options = optimset; 
my_options.Display = 'iter'; 
my_options.TolFun = 1e-3; 
my_options.MaxIter = 500; 

[x fval exitflag output] = fminunc(f, x0, my_options); 

%[x,fval] = fminunc(f,x0);
%[ rmse ] = calcVasicekNSquareRMSE( x0, cmt);