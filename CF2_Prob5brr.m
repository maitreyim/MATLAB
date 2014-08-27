%PROBLEM 5: SIMULATION OF GEOMETRIC BROWNIAN MOTION
% PARAMETERS OF r (Interest Rate),Sigma(Volatility) and
%S0(Initial Stock Price) are given in the problem as bellow
r = 0.04; 
sigma = 0.2; 
S0 = 88;    
N = 1000; % no of iterations/normal numbers
n= 10;    % no of divisions 
%Part A, Problem 5: GENERATE stock price at n time periods
for i=1:n+1
    z = randn(N,1);               % normal random number generation
    W = sqrt(i-1)*z; %generate brownian motion
    S = S0*exp(sigma*W + (r-(sigma^2)/2)*(i-1)); %generate stock price
    E(i) = mean(mean(S)); %Expected Value or Mean, This is E(Sn)
end
%plot graph
t = 1:100:N+1;
plot(t, E);
Title('Stock price Dynamics with Time');
Xlabel('Time');
YLabel('Stock Price Dynamics');
hold all %hold them all for the time being 
%since we want to insert graphs from part B as well
%PART B, Problem 5: Simulation of 6 Paths
p = 6; % no of paths(p) to be generated as asked
T = 10;     % final time as given
S(1) = S0;  % Initial Stock Price As Given
m=1:1:N+1; % x - axis of graph
% We need to generate Stock prices at 1/N time diff from each other
%Lets start with no of paths
for i=1:p
    z= randn(N,1); % Generation of Normal Random 
    %generate Stock prices at 1/N time diff till final Time
    for k=1:N
        W(k) = sqrt(1/N)*z(k); %generate borownian motion
        S(k+1) = S(k)*exp(sigma*W(k) + (r-(sigma^2)/2)*1/N);  % stock price
    end
    a1=mean(mean(S)); % take mean
    % plot S or STock Dynamics
    figure(1);
    plot(m,S);
    hold all;%hold them all in one graph
end
hold off % plot them all NOW!!
% %Question d
% r = 0.04; 
% sigma = 0.3; 
% S0 = 88;    
% N = 1000; % no of iterations/normal numbers
% n= 10;    % no of divisions in time period
% % %Part A: GENERATE stock price at n time periods
% % for i=1:n+1
% %     z = randn(N,1);               % normal random number generation
% %     W = sqrt(i-1)*z; %generate brownian motion
% %     S = S0*exp(sigma*W + (r-(sigma^2)/2)*(i-1)); %generate stock price
% %     E(i) = mean(mean(S)); %Expected Value or Mean, This is E(Sn)
% % end
% % %plot graph
% % t = 1:100:N+1;
% % plot(t, ES);
% % Title('Future Stock price and Time');
% % Xlabel('Time');
% % YLabel('Stock Price');
% % hold all
% % %PART B: Simulation of 6 Paths
% % paths = 6; % no of paths to be generated
% % T = 10;     % final time
% % S(1) = S0;  % Initial Stock Price As Given
% % m=1:1:N+1; % for plotting graph x - axis
% % %figure(2);
% % % e need to generate Stock prices at 1/N time diff from each other
% % %Lets start with no of paths
% % for i=1:paths
% %     z= randn(N,1); % Generation of Normal Random 
% %     %generate Stock prices at 1/N time diff till final Time
% %     for j=1:N
% %         W(j) = sqrt(1/N)*z(j); %generate borownian motion
% %         S(j+1) = S(j)*exp(sigma*W(j) + (r-(si^2)/2)*1/N);  % stock price
% %     end
% %     a1=mean(mean(S)); % take mean
% %     % plot graph
% %     figure(2);
% %     plot(m,S);
% %     hold all;
% % end
% % hold off