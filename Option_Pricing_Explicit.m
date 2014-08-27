function  V = Option_Pricing_Explicit(Vol, Rate, PType, Strike, Maturity, OptType, NAS)
 
% PType:   Call 'C' or Put 'P'
% OptType: European 'E' or American 'A'
% NAS:     the number of asset steps

% some buffers
S      = zeros(NAS+1, 1);          % Asset Array
Vold   = zeros(NAS+1, 1);         
Vnew   = zeros(NAS+1, 1);
Payoff = zeros(NAS+1, 1);          % payoff array
V      = zeros(NAS+1, 6);          % store option values and Greeks at current time

dS  = 2 * Strike / NAS;            % upper boundary of the stock price is twice the strike 
dt  = 0.9 / Vol ^ 2 / NAS ^ 2;     % for stability
NTS = floor(Maturity / dt) + 1;    % number of time steps
dt  = Maturity / NTS;              % to ensure that maturity is an integer number of time steps

q   = 1;
if PType == 'P'    % test for call or put
    q = -1;
end

% Terminal condition
for i = 0:NAS
    S(i+1)      = i * dS; 
    Vold(i+1)   = max(q * (S(i+1) - Strike), 0); % set up payoff
    Payoff(i+1) = Vold(i+1);
    V(i+1, 1)   = S(i+1);         % first column is the stock price
    V(i+1, 2)   = Payoff(i+1);    % second column is the payoff
end

% time loop
for k = 2:NTS+1       % time loop
    for i = 2:NAS     % asset loop. End points treated separately
        Delta = (Vold(i+1) - Vold(i-1)) /(2*dS);   % central difference
        Gamma = (Vold(i+1) - 2*Vold(i) + Vold(i-1))/(dS)^2;
        Theta = -0.5*Vol^2*S(i)^2*Gamma - Rate*S(i)*Delta + Rate*Vold(i); % Black-Scholes
        
        Vnew(i) = Vold(i) - dt*Theta;
    end
    Vnew(1)   = Vold(1)*(1 - Rate);  % Boundary condition at S = 0
    Vnew(NAS+1) = 2*Vnew(NAS) - Vnew(NAS-1); % Boundary condition at S = infinity
    
    for i = 1:NAS+1    % replace old with new
        Vold(i) = Vnew(i);
    end
    
    if OptType == 'A'   % American option, check for early exercise
        for i = 1:NAS+1
            Vold(i) = max(Vold(i), Payoff(i));
        end
    end
end

for i = 2:NAS
    V(i, 3) = Vold(i); % third column is option values
    V(i, 4) = (Vold(i+1) - Vold(i-1)) /(2*dS); % delta
    V(i, 5) = (Vold(i+1) - 2*Vold(i) + Vold(i-1))/(dS)^2; % Gamma
    V(i, 6) = -0.5*Vol^2*S(i)^2*V(i, 5) - Rate*S(i)*V(i, 4) + Rate*Vold(i); % theta
end
V(NAS+1, 3) = Vold(NAS+1);
V(1, 4)     = (Vold(2) - Vold(1)) /(dS); 
V(NAS+1, 4) = (Vold(NAS+1) - Vold(NAS)) /(dS);
V(1, 6)     = Rate*Vold(1);
V(NAS+1, 6) = -0.5*Vol^2*S(NAS+1)^2*V(NAS+1, 5) - Rate*S(NAS+1)*V(NAS+1, 4) + Rate*Vold(NAS+1);


