function [] = CF4_6A(AssetP, Strike, RiskFree, Time, Vol, nSteps)
dt = Time / nSteps;                     % Allocates the time steps
% Specifies the cost of carry (r - D)
%RR = exp(RiskFree * dt);
Up = exp(Vol * sqrt(2 * dt));           % The magnitude of an up movement
Down = 1 / Up;                          % The magnitude of a down movement
%%% Specifies the probability of up, down and mid moves for trinomial tree
P_up = ((exp(dt / 2) - exp(-Vol * sqrt(dt / 2))) / (exp(Vol * sqrt(dt / 2)) - exp(-Vol * sqrt(dt / 2)))) ^ 2;
P_down = ((exp(Vol * sqrt(dt / 2)) - exp(dt / 2)) / (exp(Vol * sqrt(dt / 2)) - exp(-Vol * sqrt(dt / 2)))) ^ 2;
P_mid = 1 - P_up - P_down;
Df = exp(-RiskFree * dt);

% u=1/d;
% d=exp((-1)*sigma*sqrt(3*delta));
% pd=(r*delta*(1-u)+(r*sigma)^2)/((u-d)*(1-d));
% pu=(r*delta*(1-d)+(r*sigma)^2)/((u-d)*(u-1));
% pm=(1-pu-pd);
% Sets up the asset movements on the trinomial tree
for i = 0:(2 * nSteps)
    State = i + 1;
    Value(State) = max(0,(AssetP * Up ^ max(i - nSteps, 0) * Down ^ max(nSteps * 2 - nSteps - i, 0) - Strike));
end

% Works backwards recursively to determine the price of the option
for tt = (nSteps - 1):-1:0
    for i = 0:(tt * 2)
        State = i + 1;
        Value(State) = (P_up * Value(State + 2) + P_mid * Value(State + 1) + P_down * Value(State)) * Df;
    end
end

Trinomial = Value(1)