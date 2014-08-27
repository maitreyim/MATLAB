% Question 3: OPTION GREEKS

epsilon = 0.02; %Is the changes in variables for different greeks
r = 0.03; % risk free interest rate as given in the question
sigma = 0.2;% Volatility as given in the question
K = 50; %Strike price as given in the question
T = 0.3846; %Time to maturity as given in the question

%Calculate the greeks for different S0 = i
S0_set = 10:2:80; % Setting Stock Price to vary from $10 to $80 with an increment of $2

delta = zeros(1, size(S0_set,2)); %Vector of deltas for different S0 
gamma = zeros(1, size(S0_set,2)); %Vector of gammas for different S0 
theta = zeros(1, size(S0_set,2)); %Vector of thetas for different S0 
vega = zeros(1, size(S0_set,2)); %Vector of vegas for different S0 
rho = zeros(1, size(S0_set,2)); %Vector of rhos for different S0 

% looping for different step sizes
for i = 1:size(S0_set,2)
    S0 = S0_set(i);

    %Delta: First Option Greek which measures Sensitivity w.r.t. an
    %incremental change in stock price
    S1 = S0; %Stock price
    S2 = S0 + epsilon; %Stock price if changed little bit
    CofS1 =option_g(S1,sigma,T,r,K ); %Price of a call for certain stock price
    CofS2 = option_g(S2,sigma,T,r,K ); %Price of a call for little change in stock price
    delta1 = (CofS2 - CofS1) / epsilon; %Calculation of delta (slope)
    delta(i) = delta1;  
    %The calculations for other greeks a very similar to the delta one.
    %Gamma: Second Option Greek which measures sensitivity of delta w.r.t
    %stock price
    S3 = S0 - epsilon;
    CofS3 = option_g(S3,sigma,T,r,K );
    gamma(i) = (CofS2 - 2*CofS1 + CofS3) / epsilon^2;    
    
%   Theta: Third Option Greek which measures the sensitivity of Call price
%   w.r.t. Change of Time
    T1 = T;
    T2 = T + epsilon;
    CofT1 = option_g(S0,sigma,T1,r,K );
    CofT2 =option_g(S0,sigma,T2,r,K );
    theta(i) = (CofT2 - CofT1) / epsilon;

    %Vega: Fourth Option Greek which measures sensitivity of Call Price
    %w.r.t. change in Volatility
    sigma1 = sigma;
    sigma2 = sigma + epsilon;
    CofSigma1 = option_g(S0,sigma1,T,r,K );
    CofSigma2 = option_g(S0,sigma2,T,r,K ); 
    vega(i) = (CofSigma2 - CofSigma1) / epsilon;

%     Rho: Fifth Option Greek which measures sensitivity of call price
%     w.r.t change in interest rate
    r1 = r;
    r2 = r + epsilon;
    CofR1 = option_g(S0,sigma,T,r1,K );
    CofR2 = option_g(S0,sigma,T,r2,K ); 
    rho(i) = (CofR2 - CofR1) / epsilon;
    
end;
% %Draw the greeks for different stock price
%subplot(2,2,1)
figure(1)
plot(S0_set , delta, '-r', S0_set , gamma, '-b', S0_set , theta, '-g', S0_set , vega, '-k', S0_set , rho, '-c', 'linewidth',2);
legend('Delta', 'Gamma', 'Theta', 'Vega', 'Rho');
Title('Option Greeks');
Xlabel('Stock Price Dynamics');
%Question 2ii:

% T_set = 0:.01:.3846;% Time Increments
% S0 = 49;% Initial Stock Price as given
% delta = zeros(1, size(T_set,2)); %Vector of deltas for different T
% for i = 1:size(T_set, 2)
%     T = T_set(i);
%     %Delta
%     S1 = S0; %Stock price
%     S2 = S0 + epsilon; %Stock price if changed little bit
%     CofS1 =option_g(S1,sigma,T,r,K ); %Price of a call for certain stock price
%     CofS2 = option_g(S2,sigma,T,r,K ); %Price of a call for little change in stock price
%     delta1 = (CofS2 - CofS1) / epsilon; %Calculation of delta (slope)
%     delta(i) = delta1;  
% end
% % Draw the change of delta w.r.t Change of Time t
% subplot(2,2,2)
% plot(T_set, delta)
% Title('Delta Dynamics w.r.t. Time');
% Xlabel('Time Dynamics');
% 
% 




