% Project 5
% get the stock values

function Price = PriceAmericanPutLSM(S0,r,si,K,T,k,poly,P,M,N,optType,S)
% Function to generate sample paths using antithetic variates
% When rnd is a vector of random numbers used to generate
% one simulated path then -rnd is used to generate another
% path.
%

Idx = zeros(M,N);
V = zeros(M,N);

% calculate essentials to make calculations easier and identifiable
dt = T/N;
disc = exp(-r*dt);
discAll = cumprod((ones(M,N).*disc),2);
flipDiscAll = fliplr(discAll);
fdiscall = flipDiscAll(1,:);

%step 2:: calc expected values at different paths and times
EV = max(K-S(:,2:end),0);

%step 3:: set Value = Expected Value at T=N
%V(:,end) = EV(:,end);
V = EV;

%step 4:: set index value to 1 if V > 0
for i=1:M
    if V(i,end)>0 
        Idx(i,end) = 1;
    end
end
if optType~='E'
    %step 5::
    for ii=N-1:-1:2
        
        % Set stock price for current column
        St=S(:,ii+1);
        % Get EV for all columns including and beyond current column
        EV_1 = EV(:,ii:end);
        % Get Index for all columns including and beyond current column
        Idx_1 = Idx(:,ii:end);
        % Get discount factors for all columns including and beyond current
        % column
        fdiscall_1 = fdiscall(1,ii:end);
        
        % Step :: 5. Estimate Y
        % multiply Expected value with Index to get intermediate Y value
        YY = EV_1.*Idx_1;
        %Get Yval
        Yval = YY*fdiscall_1';
        %Create index
        
        A = ones(k, k);
        for mm=1:k
            %build Polynomial Matrices
            f(:,mm) = P{mm}(St);
        end
        
        % build matrix A
        for j=1:k
            for i=1:k
                %             disp(f(:,j)');
                %             disp(f(:,i));
                A(i,j) = f(:,j)'*f(:,i);
                %             disp(A);
            end
            %build matrix Y
            Y(j) = Yval'*f(:,j);
        end
        
        % Form the design matrix
        %     X1 = [ones(size(St'))  St'/St' ];
        %     X2 = [ones(size(St'))  St'/St'  St'];
        %
        %     a1 = X1\Y;
        %     a2 = X2\Y;
        % get values of a 's
        a = A\Y';
        
        % find all continuation values based on a's
        %     disp(a);
        CV = zeros(M,1);
        for i=1:k
            CV = CV+ a(i)*(P{i}(St));
        end
        %     CV1 = [ones(size(St'))  St'/St' ]*a1;
        %     CV2 = [ones(size(St'))  St'/St' St']*a2;
        % create index based on Expected value being greater than Continuation
        % Value
        Idx(:,ii) = EV(:,ii) > CV;
        V(:,ii) = max(EV(:,ii), CV);
        for i=1:M
            if (Idx(i,ii)==1)
                Idx(i,ii+1:end) = 0;
            end
        end
    end
end
% Calculate discounted value of Expected Value at each path
DCF = (V.*Idx)*discAll(1,:)';
%Calculate the price of Option
Price = mean(DCF);
end