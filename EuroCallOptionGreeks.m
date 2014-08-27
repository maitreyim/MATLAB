
%                                                       r, % rate of interest
%                                                       si, % standard deviation
%                                                       T,  % time to expiry
%                                                       K   % strike price 
%                                                       S0strt start S0
%                                                       S0end end S0

function [ delta, gamma, theta, vega, rho ] = EuroCallOptionGreeks( r, sigma, T,K,S0strt, S0end)  
%EuroCallOptionGreeks Compute the greeks for European Call option prices

S0=S0strt:S0end;

% define d1 and d2
d1 = (log(S0/K) + (r+(sigma.^2))*T/2)/sigma*sqrt(T);
d2 = d1 - sigma*sqrt(T);

% calculate greeks based on formula given in lecture notes
delta = getN(d1); 
gamma = ((S0*sigma.*sqrt(T)).^(-1))*(sqrt(2*pi)).^(-1).*exp(-1*(d1.^2)/2);
theta = -S0*sigma*(sqrt(2*pi)).^(-1).*exp(-1*(d1.^2)/2).*(2*sqrt(T))^(-1) - r*K.*exp(-r*T)*N(d2);
vega = S0*sqrt(T)*(sqrt(2*pi)).^(-1).*exp(-1*(d1.^2)/2);
rho = K*T*exp(-r*T)*N(d2);

%plot graph
figure(1);
plot(S0, delta);
title('Greeks vs Stock Price');
xlabel('Initial Stock Price --> ');
ylabel('Greeks -->');
hold all
plot(S0, gamma);
hold all;
plot(S0, theta);
hold all;
plot(S0, vega);
hold all;
plot(S0, rho);
legend ('delta', 'gamma', 'theta', 'vega', 'rho');
hold off;
end

