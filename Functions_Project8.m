function CFProject8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% COMPUTATIONAL FINANCE PROJECT 8
%PRICING FIXED INCOME SECURITIES USING VASICEK/CIR/G2++ MODELS
%SUBMITTED BY MAITREYI MANDAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETER INITIALIZATION FOR QUESTION 1
nominal = 1000;
coupon = 30;
r0 = 0.05; %interest rate at t = 0;
sigma = 0.1; %Volatility
kappa = 0.82; %Speed of reversion
aveR = 0.05; %average r;
S = 0.5; % Time to maturity of the bond, S > T
M = 1000; %How many r paths to simulate for bond prices 
T = 3/12; %Time to maturity of the options
no = ceil(365*T); %How many time steps in each interest rate path(daily)
Mo = 100; %How many r paths to simulate for option prices
K = 950; %Strike price
oPrices = zeros(100,1); %option prices for each Mo path
%%%%%%%%%%%Question 1a,1b,1c,1d,1e%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%PART 1 A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the price of a pure discount bond using Vasicek
price = PriceDiscountBond(nominal, 0, 0, r0, 0, 0, sigma, 0, kappa,.... 
aveR, 0, 0, S, M, 'Vasicek')
%%%%%%%%%%%%%PART 1 B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Get the price of a coupon paying bond using Vasicek
% price = PriceCouponBond(1000, 30, r0, sigma, kappa, aveR, 0, 4, M,....... 
% 0.5, 'Vasicek') 
% % 
% %%%%%%%%%%%PART1C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Number of paths to get the value of the option
% % for i = 1:1000
% % % %     %Create each path using Vasicek model
% % ro = OneFactorVasicek(r0, sigma, kappa, aveR, T, no, 1);
% % % %     %Get the price of a discount bond for each path
% P = PriceDiscountBond(nominal, 0, 0, ro(no), 0, 0, sigma, 0,............. 
% kappa, aveR, 0, 0, S-T, M, 'Vasicek');
% % % %     %Get the option price for each path
% % oPrices(i) = exp(-mean(ro)*T)*max(P-K,0);
% % end;
% % option = mean(oPrices)
% %Get the average option price
% % 
% %%%%%%%%%%%PART1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Number of paths to get the value of the option
% % for i = 1:1000
% % % %     %Create each path using Vasicek model
% % ro = OneFactorVasicek(r0, sigma, kappa, aveR, T, no, 1);
% % % %     %Get the price of a coupon paying bond for each path
% P = PriceCouponBond(nominal, coupon, ro(no), sigma, kappa, aveR, T,...... 
% 4.0, M, 0.5, 'Vasicek');
% % % %     %Get the option price for each path
% % oPrices(i) = exp(-mean(ro)*T)*max(P-K,0);
% % end;
% % option = mean(oPrices)%Get the average option price
% % 
% %%%%%%%%%%PART 1E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Number of paths to get the value of the option
% % for i = 1:100
% % % % %     %Create each path using Vasicek model
% % ro = OneFactorVasicek(r0, sigma, kappa, aveR, T, no, 1);
% % % %Get the price of a coupon paying bond for each path using the explicit method
% P = ExplicitVasicekPriceCouponBond(nominal, coupon, ro(no), sigma,....... 
% kappa, aveR, T, 4.0, 0.5);
% % % % %     %Get the option price for each path
% %  oPrices(i) = exp(-mean(ro)*T)*max(P-K,0);
% % end;
% % option = mean(oPrices);%Get the average option price
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%Question 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%QUESTION 2 A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Same input parameters as in question 1, just different values
% % r0 = 0.05; sigma = 0.12; kappa = 0.92; aveR = 0.055;
% % S = 1; T = 0.5; no = ceil(T*365); Mo = 100; K = 950;
% % nominal = 1000; M = 100;
% % % % %Number of paths to get the value of the option
% % for i = 1:10000
% % % %     %Create each path using CIR model
% % ro = OneFactorCIR(r0, sigma, kappa, aveR, T, no, 1);
% % % %     %Get the price of a pure discount bond for each path
% P = PriceDiscountBond(nominal, 0, 0, ro(no), 0, 0, sigma, 0, kappa,...... 
% aveR, 0, 0, S-T, M, 'CIR');
% % % %     %Get the option price for each path
% % oPrices(i) = exp(-mean(ro)*T)*max(P-K,0);
% % end;
% % option = mean(oPrices)%Get the average option price
% % % 
% %%%%%%%%%%%QUESTION 2C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %T is the maturity of the option, S is the maturity of the bond.
% % T = 0.5; S = 1; t = 0; %t is the initial time in the calculation
% % % % %Get the explicit price
% % price = ExplicitCIRPriceOption(1, 0.950, 0.05, 0.12, 0.92, 0.055, T, S, t)
% % 
% %%%%%%%%%%%Question 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%PARAMETER INITIALIZATION FOR G2++ model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 = 0; y0 = 0; phi0 = 0.03; r0 = 0.03; sigma = 0.03; nu = 0.08; a = 0.1;...
%     b = 0.3;
% rho = 0.7; S = 1; T = 0.5; no = ceil(T*365); Mo = 100; K = 900;
% nominal = 1000; M = 1000;
% %%Number of paths to get the value of the option
% for i = 1:1000
% % %     %Create each path using G2++ model
% ro = TwoFactorG2(x0, y0, r0, phi0, rho, sigma, nu, a, b, T, no, 1);
% % %     %Get the price of a discount bond for each path
% P = PriceDiscountBond(nominal, x0, y0, ro(Mo), phi0, rho, sigma, nu, 0, 0,... 
%     a, b, S-T, M, 'G2++');
% % %     %Get the option price for each path, this is put.
% oPrices(i) = exp(-mean(ro)*T)*max(K-P,0);
% end;
% option = mean(oPrices);%Get the average option price
% % 
% % %Get the explicit price of the option price using G2++
% optionPrice = ExplicitG2PriceOption('P', 0.9, rho, sigma, nu, phi0, x0,.... 
%     y0, a, b, 0, T, S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
end
function [ curve ] = OneFactorCIR( r0, sigma, kappa, aveR, T, n, M )
%Function to generate stochastic interest rate using CIR model
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
r = zeros(1,M);%Average rate of each path

for j = 1:M%To create each path
    Z = randn(n, 1);%Create random numbers for each path
    for i = 1:n-1%To create each point in the path
        deltaR = kappa*(aveR - rt(i))*deltaT + sigma*sqrt(rt(i))*sqrt(deltaT)*Z(i,1);
        rt(i+1) = rt(i) + deltaR;
    end;
    if(M ~= 1)%If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;
curve = r;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ curve ] = OneFactorVasicek( r0, sigma, kappa, aveR, T, n, M )
%Function to generate stochastic interest rate using Vasicek model
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
deltaT = T/n;%Size of each time step
rt = zeros(n,1);%Vector to capture all the points in one path of interest rate
rt(1) = r0;%The first point in the path
r = zeros(1,M);%Average rate of each path

for j = 1:M %To create each path
    Z = randn(n,1); %Create random numbers for each path
    for i = 1:n-1 %To create each point in the path
        deltaR = kappa*(aveR - rt(i))*deltaT + sigma*sqrt(deltaT)*Z(i,1);
        rt(i+1) = rt(i) + deltaR;
    end;
    if(M ~= 1) %If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;

curve = r;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = PriceDiscountBond(nominal, x0, y0, r0, phi0, rho, sigma, nu, kappa, aveR, a, b, T, M, model)
%Function to calculate the price of a discount bond with stochastic
%interest rate. To get the evolution of the interest rate we can use three
%different models. Vasicek, CIR or G2++
%Nominal is the face value of the bond.
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths.
%Parameters only for G2++ model: x0, y0, phi0, rho, nu, a and b

n = ceil(365*T);%How many time steps(days)
%What interest rates model to use
if(strcmp(model,'CIR'))
    r = OneFactorCIR(r0, sigma, kappa, aveR, T, n, M);
elseif(strcmp(model,'G2++'))
    r = TwoFactorG2(x0, y0, r0, phi0, rho, sigma, nu, a, b, T, n, M);
else
    r = OneFactorVasicek(r0, sigma, kappa, aveR, T, n, M);
end;
price = nominal*mean(exp(-r*T)*1);%Caluclate the price of a discount bond
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = PriceCouponBond(nominal, coupon, r0, sigma, kappa, aveR, Tstart, Tend, M, compounded, model)
%Function to calculate the price of a coupon paying bond with stochastic
%interest rate. To get the evolution of the interest rate we use the Vasicek model
%Nominal is the face value of the bond. Coupon is how much the coupon
%payment is.
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%Tstart is the time we want to get the price at.
%Tend is the maturity of the bond, n is number of time step in one path
%M is number of paths, compounded is how often we compound each year.

PV = zeros(8,1);%Present value of each payment
reminder = rem(Tend-Tstart,compounded);
count = 1;
from = Tstart;%Time we want to get the value of the bond
to = Tend - reminder;%From start time, how long time it is to maturity of the bond

for Ti = from:compounded:to
   if(Ti == to)
       %Get the present value of the final payment
       PV(count) = PriceDiscountBond(nominal+coupon, 0, 0, r0, 0, 0, sigma, 0, kappa, aveR, 0, 0, Ti, M, model);
   elseif(Ti == 0)
       count = 0;
   else
       %Get the present value of the coupon payments
       PV(count) = PriceDiscountBond(coupon, 0, 0, r0, 0, 0, sigma, 0, kappa, aveR, 0, 0, Ti, M, model);
   end;
   count = count + 1;
end;
%Sum up all the payments
price = sum(PV);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ curve ] = TwoFactorG2( x0, y0, r0, phi0, rho, sigma, nu, a, b, T, n, M )
%Function to generate stochastic interest rate using G2++ model
%This is a two factor model where the interest rates are generated from two
%factors.
%x0 and y0 are the initial values of the two factors
%phiO is the initial value of the shift funtion phit.
%r0 is the initial interest rate, sigma and nu are the volatilities for xt and yt,
%a and b are the drifts for xt and yt
%rho is the correlation between the two factors
%T is the maturity of the bond, n is number of time step in one path
%M is number of paths. n is the number of steps.

deltaT = T/n;%Size of each time step
%Vectors to capture all the points in one path of interest rate
rt = zeros(n,1); xt = zeros(n,1); yt = zeros(n,1); phit = zeros(n,1);
%The first point in the path
rt(1) = r0; xt(1) = x0; yt(1) = y0; phit(1) = phi0;
r = zeros(1,M);%Average rate of each path

for j = 1:M%To create each path
    Z = randn(n,2);%Create random numbers for each path, both xt and yt
    for i = 1:n-1%To create each point in the path
        %We need to find the changes in xt, yt and phit
        dx = -a*xt(i)*deltaT+sigma*sqrt(deltaT)*Z(i,1);
        dy = -b*yt(i)*deltaT+nu*sqrt(deltaT)*(rho*Z(i,1)+sqrt(1-rho^2)*Z(i,2));
        dphi = 0;
        %Get the next point for xt, yt and phit
        xt(i+1) = xt(i) + dx;
        yt(i+1) = yt(i) + dy;
        phit(i+1) = phit(i) + dphi;%phi is a shift function
        %Get each interest rate point, r, from xt, yt and phit
        rt(i+1) = xt(i+1) + yt(i+1) + phit(i+1);
    end;
    if(M ~= 1)%If we need to use more than one path we get the average interest rate for each path
        r(j) = mean(rt);
    else %Else we return the hole path.
        r = rt;
    end;
end;

curve = r;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = ExplicitCIRPriceDiscountBond(nominal, rt, sigma, kappa, aveR, T, t)
%Function to calculate the price of a discount bond with stochastic
%interest rate using explicit CIR formula. 
%Nominal is the face value of the bond.
%rt is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path
%t is the initial time

%Necessary parameters to calculate the explicit price
h1 = sqrt(kappa^2+2*sigma^2);
h2 = (kappa+h1)/2;
h3 = 2*kappa*aveR/sigma^2;
B = (exp(h1*(T-t))-1)/(h2*(exp(h1*(T-t))-1)+h1);
A = ((h1*exp(h2*(T-t)))/(h2*(exp(h1*(T-t))-1)+h1))^h3;
%Price of a pure discount bond using CIR explicit formula
price = nominal*A*exp(-B*rt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = ExplicitCIRPriceOption(nominal, K, rt, sigma, kappa, aveR, T, S, t)
%Function to calculate the price of a European Call option on pure discount bond 
%with stochastic interest rate using explicit CIR formula. 
%Nominal is the face value of the bond. K is the strike price
%sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the option, S is the maturity of the bond
%t is the initial time from where we want to get the price of the option

theta = sqrt(kappa^2+2*sigma^2);
phi = 2*theta/(sigma^2*(exp(theta*(T-t))-1));
chi = (kappa+theta)/sigma^2;

%To get A and B for S-T
h1 = sqrt(kappa^2+2*sigma^2);
h2 = (kappa+h1)/2;
h3 = 2*kappa*aveR/sigma^2;
B = (exp(h1*(S-T))-1)/(h2*(exp(h1*(S-T))-1)+h1);
A = ((h1*exp(h2*(S-T)))/(h2*(exp(h1*(S-T))-1)+h1))^h3;

rStar = log(A/K)/B;
%Parameters to by used in the Chi-Squared CDF.
x1 = 2*rStar*(phi+chi+B);
p1 = (4*kappa*aveR)/sigma^2;
q1 = (2*phi^2*rt*exp(theta*(T-t)))/(phi+chi+B);

x2 = 2*rStar*(phi+chi);
p2 = (4*kappa*aveR)/sigma^2;
q2 = (2*phi^2*rt*exp(theta*(T-t)))/(phi+chi);

%Get the value of the bonds with maturity at S and T
PS = ExplicitCIRPriceDiscountBond(nominal, rt, sigma, kappa, aveR, S, 0);
PT = ExplicitCIRPriceDiscountBond(nominal, rt, sigma, kappa, aveR, T, 0);

%Explicit price of the option
price = 1000*(PS*ncx2cdf(x1,p1,q1)-K*PT*ncx2cdf(x2,p2,q2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = ExplicitG2PriceDiscountBond(sigma, nu, rho, phi, a, b, xt, yt, T, t)
%Function to calculate the explicit pure discount bond price using G2++ model
%This is a two factor model where the interest rates are generated from two
%factors.
%xt and yt are two factors
%phit is shift funtion.
%sigma and nu are the volatilities for xt and yt,
%a and b are the drifts for xt and yt
%rho is the correlation between the two factors
%T is the maturity of the bond, t is the initial time

%To calculate V I split it into three parts x,y and z
x = (sigma^2/a^2)*((T-t)+(2*exp(-a*(T-t))/a)-(exp(-2*a*(T-t))/(2*a))-(3/(2*a)));
y = nu^2/b^2*(T-t+2/b*exp(-b*(T-t))-1/(2*b)*exp(-2*b*(T-t))-3/(2*b));
z = 2*rho*sigma*nu/(a*b)*(T-t+(exp(-a*(T-t))-1)/a+(exp(-b*(T-t))-1)/b-(exp(-(a+b)*(T-t))-1)/(a+b));
%Get V
V = x + y + z;

%Also splitted the price calculation into three parts f, s and t
f = phi*(T-t);
s = (1-exp(-a*(T-t)))*xt/a;
t = (1-exp(-b*(T-t)))*yt/b;
%Get the price
price = exp(-f-s-t+V/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = ExplicitG2PriceOption(kind, K, rho, sigma, nu, phi, xt, yt, a, b, t, T, S)
%Function to calculate the explicit option price using G2++ model
%This is a two factor model where the interest rates are generated from two
%factors.
%Kind is Call (C) or Put (P)
%xt and yt are two factors
%phit is shift funtion.
%sigma and nu are the volatilities for xt and yt,
%a and b are the drifts for xt and yt
%rho is the correlation between the two factors
%T is the maturity of the option, S of the bond and t is the initial time

%Three parts to calculate SIGMA^2
x = sigma^2/(2*a^3)*(1-exp(-a*(S-T)))^2*(1-exp(-2*a*(T-t)));
y = nu^2/(2*b^3)*(1-exp(-b*(S-T)))^2*(1-exp(-2*b*(T-t)));
z = 2*rho*sigma*nu/(a*b*(a+b))*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))*(1-exp(-(a+b)*(T-t)));
SIG2 = x + y + z;
%Get two bond prices with different maturity days
PS = ExplicitG2PriceDiscountBond(sigma, nu, rho, phi, a, b, xt, yt, S, t);
PT = ExplicitG2PriceDiscountBond(sigma, nu, rho, phi, a, b, xt, yt, T, t);
%Call or Put
if(kind == 'C')
    X1 = log(PS/(K*PT))/sqrt(SIG2)+sqrt(SIG2)/2;
    X2 = log(PS/(K*PT))/sqrt(SIG2)-sqrt(SIG2)/2;
    oPrice = PS*normcdf(X1,0,1)-PT*K*normcdf(X2,0,1);
elseif(kind == 'P')
    X1 = log(K*PT/PS)/sqrt(SIG2)-sqrt(SIG2)/2;
    X2 = log(K*PT/PS)/sqrt(SIG2)+sqrt(SIG2)/2;
    oPrice = -PS*normcdf(X1,0,1)+PT*K*normcdf(X2,0,1);
else
    oPrice = 0;
end;
%Get the option price
price = 1000*oPrice;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ price ] = ExplicitVasicekPriceCouponBond(nominal, coupon, r0, sigma, kappa, aveR, Tstart, Tend, compounded)
%Function to calculate the price of a coupon paying bond with stochastic
%interest rate using explicit Vasicek formula.
%Nominal is the face value of the bond. Coupon is how much the coupon
%payment is.
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%Tstart is the time we want to get the price at.
%Tend is the maturity of the bond, n is number of time step in one path

PV = zeros(8,1);%Present value of each payment
reminder = rem(Tend-Tstart,compounded);
count = 1;
from = Tstart;%Time we want to get the value of the bond
to = Tend - reminder;%From start time, how long time it is to maturity of the bond

for Ti = from:compounded:to
   if(Ti == to)
       %Get the present value of the final payment
       PV(count) = ExplicitVasicekPriceDiscountBond(nominal+coupon, r0, sigma, kappa, aveR, Ti);
   elseif(Ti == 0)
       count = 0;
   else
       %Get the present value of the coupon payments
       PV(count) = ExplicitVasicekPriceDiscountBond(coupon, r0, sigma, kappa, aveR, Ti);
   end;
   count = count + 1;
end;
%Sum up all the payments
price = sum(PV);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ price ] = ExplicitVasicekPriceDiscountBond(nominal, r0, sigma, kappa, aveR, T)
%Function to calculate the price of a discount bond with stochastic
%interest rate using explicit Vasicek formula. 
%Nominal is the face value of the bond.
%r0 is the initial interest rate, sigma is the volatility,
%kappa is the speed of reversion, aveR is the long term mean
%T is the maturity of the bond, n is number of time step in one path

B = 1/kappa*(1-exp(-kappa*T));
A = exp((aveR-sigma^2/(2*kappa^2))*(B-T)-sigma^2/(4*kappa)*B^2);
%Price of the pure discount bond
price = nominal*A*exp(-B*r0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





