function [P] = AmCallImplicit_2(OptionType,SO,K,r,q,T,sig,Smin,Smax,Ds,Dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function prices a Vanilla European Call/Put using the Implicit Scheme of   %                  
% the Finite Difference Method.                                                   %           
% Parameters are as follows:                                                      %          
% OptionType = 1 for a Call or 0 for a Put                                        %              
% SO = initial asset price ; K = strike price ; r = risk free rate ; q = dividend %
% rate ; T = time to maturity ; sig = volatility ; Smin = minimum stock price ;   %   
% Smax = maximum stock price ; Ds = stock price step size ; Dt = time step        %           
% size.                                                                           %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate number of stock price steps and take care of rounding.
N = round((Smax - Smin) / Ds);
Ds = (Smax - Smin) / N;
% Calculate number of time steps and take care of rounding.
M = round(T/Dt);
Dt = (T/M);

MI=zeros(N,N); % MI matrix
S=zeros(N,1); % stock price vector
V=zeros(N,1); % option value vecto
matsol=zeros(N,M+1); % solution matrix

for i=1:1:N % Generate S and V vectors
    S(i)=Smin + i*Ds;
    if OptionType == 1
        V(i)=max(S(i)-K,0); % Call: Payoff that is initial condition
    else
        V(i)=max(K-S(i),0); % Put: Payoff that is initial condition
    end    
end

for i=1:1:N % Build MI matrix
    % Set up coefficients
    Alpha = 0.5*(sig^2)*(S(i)^2)*(Dt/(Ds^2));
    Beta = (r-q)*S(i)*(Dt/(2*Ds));
    Bdi=-Alpha+Beta;
    Di=1+r*Dt+2*Alpha;
    Adi=-Alpha-Beta;
    % Fill MI matrix
    if i==1
        MI(i,i) = Di + 2*Bdi;
        MI(i,i+1) = Adi - Bdi;
    elseif i==N
        MI(i,i-1) = Bdi - Adi;
        MI(i,i) = Di + 2*Adi;
    else
        MI(i,i-1) = Bdi;
        MI(i,i) = Di;
        MI(i,i+1) = Adi;
    end
end

matsol(:,1)=V; % Initiate first column of matrix solution with payoff that
               % is initial condition

invMI = MI^-1; % Invert matrix MI before performing calculations

for k=1:M % Generate solution matrix
    matsol(:,k+1)=invMI*matsol(:,k);
end

% find closest point on the grid and return price
% with a linear interpolation if necessary

DS = SO-Smin;

indexdown = floor(DS/Ds);
indexup = ceil(DS/Ds);

if indexdown == indexup
    P = matsol(indexdown,M+1);
else
    P = matsol(indexdown,M+1)+ (SO - S(indexdown))/(S(indexup)-S(indexdown))...
       *(matsol(indexup,M+1) - matsol(indexdown,M+1));
end

