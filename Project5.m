% Project 5
% *********Q 1 a b c d************
% get the stock values
clear all; close all; clc;
S0 = [36 40 44]; si = 0.2; r=0.06; 
%M=100000; N = 100;
M=1000;N=100;
K = 40; T = [0.5 1 2]; k = [2 3 4];
p = {'Hermite'; 'Laguerre'; 'Monomial'};

% get the basis functions for least squares estimation
% get Laguerre
H = Hermite(); L = Laguerre(); MO = Monomial();
H = H(1:5); L = L(1:5); MO = MO(1:5);
P(:,1) = H; P(:,2) = L; P(:,3) = MO; 

% Binomial American Put Price function from Project 4 for comparison 
% set parameters
treeTyp='B'; CallOrPut='P'; logPrice='N'; % type of options
optTyp='A'; % calc American

% call function to get American Put price
% need four loops to get a combination of all prices:
ctrP=1;
Price{ctrP,1}={'Initial Stock Price'}; Price{ctrP,2}={'Time to Expiry'};
Price{ctrP,3}={'Type of PolyNomial'};
Price{ctrP,4}= {'No of polynomials'}; Price{ctrP,5}={'Price LSMC'};
Price{ctrP,6}= {'American Put Price Binomial'}; 
%Price{ctrP,6}={'BS Put Price'};
for i=1:3   %loop for Stock
    for j=1:3   %loop for Time
        % calculate price for comparison with standard functions
        % [callBS putBS]= blsprice(S0(i),K,r,T(j),si);
        [u,d,prob] = Q1b_bino(si,r,T(j)/N); % get u d p
        [cpPrice] = compBiTriOptionPr( r, si,S0(i),K, T(j), u, d, ...
            prob, 1-prob,0, N, treeTyp, optTyp, CallOrPut,logPrice);
        %step 1:: Calc Stock prices at different paths and times
        S = StockPathsAntithetic(S0(i),r,si,T(j)/N,N,M)';
        for l=1:3   %loop for k
            for m=1:3 %loop for polynomial to use
                ctrP = ctrP+1;
                Price{ctrP,1}=S0(i);Price{ctrP,2}=T(j);
                Price{ctrP,3}=p{m};Price{ctrP,4}=k(l);
                Price{ctrP,5} = PriceAmericanPutLSM(S0(i),r,si,K,T(j),k(l),p{m},P(:,m),M,N,'A',S);
                % Price{ctrP,6} = putBS; 
                Price{ctrP,6} = cpPrice;
                disp(['American Put Price Using LSM = ' num2str(Price{ctrP,5}) ...
                    ' :: Using S0:' num2str(S0(i)) ' :: Time:' num2str(T(j)) ...
                    ' :: k:' num2str(k(l)) ' Polynomial:' p{m} ...
                    ' :: Num of Paths:' num2str(M) ' :: Time Steps:' num2str(N)]);
                
            end
            disp('--------------------------------------------');
        end
        disp('--------------------------------------------');
    end
    disp('--------------------------------------------');
end
% S = StockPathsAntithetic(S0(i),r,si,T(1)/N,N,M)';
% Price = PriceAmericanPutLSM(S0(1),r,si,K,T(1),k(1),p(2),P(:,2),M,N,'A',S);
%[callBsPrice putBsPrice]= blsprice(S0(1),K,r,T(2),si);

% % ****************Q 2 *****************
% clear all; close all; clc;
% S0 = [36 40 44]; si = 0.2; r=0.06; 
% K = 40; T = [0.5 1 2];
% %M=100000;N=100;
% M=100000;N=100;
% 
% ctr=1;
% Q2Price{ctr,1}={'Initial Stock Price'}; Q2Price{ctr,2}={'Time to Expiry'};
% Q2Price{ctr,3}={'European Put Option Price : MC'}; 
% Q2Price{ctr,4}={'Black Scholes Put Option Price'}; 
% for i=1:3 % initial stock price loop
%     for j=1:3 % Time to expiry loop
%         ctr=ctr+1;
%         Q2Price{ctr,1}=S0(i);Q2Price{ctr,2}=T(j);
%         Q2Price{ctr,3} = PriceAmericanPutLSM(S0(i),r,si,K,T(j),1,'Laguerre',1,M,N,'E');
%         [callBSPrice, Q2Price{ctr,4}]=blsprice(S0(i),K,r,T(j),si);
%     end
% end

% %*****************Q 3 a b*****************
% clear all; close all; clc;
% S0 = 65; si = 0.2; r=0.06; t=0.2; T=1;
% %M=100000;N=100;
% M=1000;N=100;
% dt=T/N;
% k = [2 3 4];
% p = {'Hermite'; 'Laguerre'; 'Monomial'};
% % get the basis functions for least squares estimation
% % get Laguerre
% H = Hermite(); L = Laguerre(); MO = Monomial();
% H = H(1:5); L = L(1:5); MO = MO(1:5);
% P(:,1) = H; P(:,2) = L; P(:,3) = MO; 
% 
% % generate Antithetic Stock Prices
% % step 1:: Calc Stock prices at different paths and times
% S = StockPathsAntithetic(S0,r,si,dt,N,M)';
% % find Strike Price
% fp = (N*t/T)+1; % forward start point
% NN = N - fp+1; % modified N = number of time steps that get reduced
% K = S(:,fp); % Strike Price
% newT = T - t; % new Time to Expiry based on Strike price calculation
% disc = exp(-r*dt*t);
% % generate Stock prices to pass to LSMC function
% SfS = S(:,fp:end);
% % use the MC Stock prices to calculate European and American put option 
% % values for different forward start strike prices calculated above
% for i=1:M
%     PriceEA(i,1) = PriceAmericanPutLSM(S0,r,si,K(i),newT,k(3),p(3),P(:,3),M,NN,'E',SfS);
%     PriceEA(i,2) = PriceAmericanPutLSM(S0,r,si,K(i),newT,k(3),p(3),P(:,3),M,NN,'A',SfS);
% end
% EuroPutOption = mean(PriceEA(:,1))*disc;
% AmerPutOption = mean(PriceEA(:,2))*disc;
% disp(['Forward-Start European Put Option Price is : ' num2str(EuroPutOption)]);
% disp(['Forward-Start American Put Option Price is : ' num2str(AmerPutOption)]);
% 












