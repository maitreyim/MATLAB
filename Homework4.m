
%Question 1:
%HW4Question1();

%Question 2:
%value = HW4Question2();

%Question 3:
%Greeks();

%Question 5;
%HW4Question5()

%Question 6:
%HW4Question6();

%Question 7:
S0 = 32; %Stock price
K = 30; %Strike price
T = 0.5; %Time to maturiy
r = 0.05; %Interest rate
sigma = 0.24; %Volatility
type = 'C';
N = 500; %Number of points
b1 = 2; %Base 1
b2 = 5; %Base 2
seqBaseB1 = GetHalton(N/2, b1); %Halton sequence with base = b1
seqBaseB2 = GetHalton(N/2, b2); %Halton sequence with base = b2
%Generate normal distributed random numbers from Halton sequences
normalSequence = BoxMuller(N, seqBaseB1, seqBaseB2);
valuesAtT = zeros(1, N); %Vector to accumulate all the stock prices at T 

%Check if all the constrains are ok
if(S0 >= 0 && K >= 0 && sigma >= 0 && T >= 0 && N >= 0)
    for t = 1:N
        %Generate N different stock prices from the random numbers above. 
        wt = sqrt(T) * normalSequence(t); 
        St = S0 * exp(sigma*wt + (r-sigma^2/2)*T); %Value of the underlying stock
        %Different calculations for Calls or Puts
        if type == 'C'
            valuesAtT(t) = max(St - K, 0);
        elseif type == 'P'
            valuesAtT(t) = max(K - St, 0);
        else
            valuesAtT(t) = 0;
        end;
    end;
    %Take the average of these stock prices to get one stock price in each point
    %and discount it to time zero.
    valueOption = exp(-r*T)*mean(valuesAtT); %Get the call value;
    [call, put] = blsprice(S0, K, r, T, sigma); %Compare to Black-Scholes value
else
    valueOption = 0; %If constrain above not ok, value = 0
end;













