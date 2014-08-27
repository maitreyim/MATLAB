%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Devdeep Sarkar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conditional Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% The Given Basic Parameters %%%%%%%%%%%%%

K = 100;                       % The strike price 
S0 = 100;
vol0 = 0.04;                    % Initial value of the stochastic vol
T = 1;
r = 0.05;
q = 0;                         % zero dividend

pv = exp(-r*T);                % Discount factor

m = 50;                        % Number of steps (N in the problem)
del = T/m;                     % Step size
n = 10000;                     % Number of simulations = 100K

%%%%%% %%%%%% %%%%% %%%%%

alpha = 0.10;                  % drift parameter for stochastic vol
psi = 0.10;                    % vol parameter of vol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Euler Scheme + Standard Monte Carlo Technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(m+1,n);
vol = zeros(m+1,n);

z1 = randn(m,n);  
z2 = randn(m,n);

S(1,:) = S0;
vol(1,:) = vol0;

for i = 1:m
    S(i+1,:) = S(i,:) + r*del*S(i,:) + sqrt(del)*(vol(i,:).^(0.5)).*S(i,:).*z1(i,:);
    vol(i+1,:) = vol(i,:) + alpha*del*vol(i,:) + psi*sqrt(del)*vol(i,:).*z2(i,:);    
end

c = pv*max(0,S(m+1,:)-K);                % Discounted payoff

optionA = mean(c)
errorA = std(c)/sqrt(n)

% Comparison with Black-Scholes:

bsPrice = blsprice(S0,K,r,T,sqrt(vol0))
diff = (bsPrice-optionA)/errorA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) Conditional Monte Carlo Technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol = zeros(m+1,n);
vol(1,:) = vol0;
z = randn(m,n); 
tot = zeros(1,n);

for i = 1:m
    vol(i+1,:) = vol(i,:) + alpha*del*vol(i,:) + psi*sqrt(del)*vol(i,:).*z2(i,:);  
    tot = tot + vol(i+1,:);
end

averageVol = tot/m;
c = blsprice(S0,K,r,T,sqrt(averageVol));

optionB = mean(c)
errorB = std(c)/sqrt(n)





