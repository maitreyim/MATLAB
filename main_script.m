r0 = 0.05; sigma = 0.12; kappa = 0.92; aveR = 0.055;
S = 1; T = 0.5; no = ceil(T*365); Mo = 100; K = 950;
nominal = 1000; M = 100;
% %Number of paths to get the value of the option
for i = 1:10000
    %Create each path using CIR model
    ro = OneFactorCIR(r0, sigma, kappa, aveR, T, no, 1);
    %Get the price of a pure discount bond for each path
    P = PriceDiscountBond(nominal, 0, 0, ro(no), 0, 0, sigma, 0, kappa, aveR, 0, 0, S-T, M, 'CIR');
    %Get the option price for each path
    oPrices(i) = exp(-mean(ro)*T)*max(P-K,0);
end;
option = mean(oPrices);%Get the average option price