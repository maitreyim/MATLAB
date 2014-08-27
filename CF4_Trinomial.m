function [Trinomial] = CF4_Trinomial(S, K, r, T, sigma, nSteps)
%Question5, PartA: European Call Pricing Using Trinomial Method
delta = T / nSteps;  % Allocates the time steps
                   
%%% Specifies the probability of up, down and mid moves for trinomial tree
d=exp((-1)*sigma*sqrt(3*delta));
u=1/d;
pd=(r*delta*(1-u)+(r*delta)^2 + (sigma^2)*delta)/((u-d)*(1-d));
pu=(r*delta*(1-d)+(r*delta)^2 + (sigma^2)*delta)/((u-d)*(u-1));
pm=(1-pu-pd);
Df = exp(-r * delta);

% Sets up the asset movements on the trinomial tree
for i = 0:(2 * nSteps)
    State = i + 1;
    Value(State) = max(0, (S * u ^ max(i - nSteps, 0) * d ^ max(nSteps * 2 - nSteps - i, 0) - K));
end

%S0 + u*uPower + d*dPower;
% Works backwards recursively to determine the price of the option
for l = (nSteps - 1):-1:0 % Going back one step
    for i = 0:(l * 2)
        State = i + 1;
        Value(State) = (pu * Value(State + 2) + pm* Value(State + 1) + pd * Value(State)) * Df;
    end
end

Trinomial = Value(1);

end
