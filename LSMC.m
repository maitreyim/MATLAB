% Pricing American Put using Longstaff-Schwartz method
% Stock underlying dynamics: dS = r*S dt + sig*S dW
clear all; close all
r = 0.06;
sig = 0.2;
T = .5;
S0 = 36;
K = 40;
N = 200;
M = 50000;
h = T/N;
% Simulate M independent paths first,for regression.
Splus = ones(M, N);
Sminus = ones(M, N);
Splus(:,1) = S0*ones(M,1);
Sminus(:,1) = S0*ones(M,1);
for n = 2:N
  dW = sqrt(h)*randn(M,1);
  Splus(:,n) = Splus(:,n-1).*exp(r*h+sig*dW);
  Sminus(:,n) = Sminus(:,n-1).*exp(r*h+sig*dW*(-1));
end
S=[Splus;Sminus];
InMat=zeros(size(S));
V=NaN (size(S));
V(:,end)=max(0,K-S(:,end));
InMat(:,end)=(V(:,end)>0);

for n=(N-1):-1:1
    itmpath=(K - S(:,n) > 0 );
    itmS=S(itmpath,n);
    X = itmS;
    regmat=[exp(-X/2), exp(-X/2).*(1-X), exp(-X/2).*( 1 - 2*X - X.^2/2), ...
        exp(-X/2).*(1 - 3*X + 3*X.^2/2 - X.^3/6)];
    % regmat=[itmS.^0,itmS.^1,itmS.^2,itmS.^3,itmS.^4];
    CV=V(itmpath,n+1)*exp(-r*h);
 
    result=regress(CV,regmat);
    estimate=regmat*result;
    % estimated pay-off 
    InMat(itmpath,n)=(estimate<(K-S(itmpath,n)));
    exercised = (InMat(:,n)==1);
    InMat(exercised,(n+1):end)=0;
    V(:,n)=exp(-r*h)*V(:,n+1);

    V(exercised,n) = (K - S(exercised, n));
    assert( all( V(:,n) >= 0 ))
end

mean(V(:, 1))

discountmat = repmat( exp(-h*r).^(0:(N-1)), 2*M, 1);
test = sum(discountmat .* (K - S) .* InMat, 2);





















