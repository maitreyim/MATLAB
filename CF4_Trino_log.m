function [Trinomial] = CF4_Trino_log(S, K, r, T, sigma, nSteps)
%Question5, PartA: European Call Pricing Using Trinomial Method
delta = (T/nSteps);  % Allocates the time steps
%%% Specifies the probability of up, down and mid moves for trinomial tree
Xu=sigma*sqrt(3*delta);
Xd= exp(-sigma*sqrt(3*delta));
pd=0.5*(((sigma^2)*delta+((r-(sigma^2)/2)^2)/(Xu^2))-((r-(sigma^2)/2)*delta)/Xd);
pu=0.5*(((sigma^2)*delta+((r-(sigma^2)/2)^2)/(Xu^2))+((r-(sigma^2)/2)*delta)/Xd);
 
pm=(1-pu-pd);
Df = (-r * delta);
% Sets up the asset movements on the trinomial tree
for i = 0:(2 * nSteps)
    State = i + 1;
   Value(State) = max(0, (S +Xu*( max(i - nSteps, 0) )+Xd *(max(nSteps * 2 - nSteps - i, 0) - K)));
end

% Works backwards recursively to determine the price of the option
for l = (nSteps - 1):-1:0 % Going back one step
    for i = 0:(l * 2)
        State = i + 1;
        
        Value(State) = (pu * Value(State + 2) + pm* Value(State + 1) + pd * Value(State)) +Df;
    end
end

Trinomial = Value(1);

end


 
    
 