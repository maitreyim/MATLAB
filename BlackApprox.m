function [ option ] = BlackApprox( r, si, S0, K, T, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BlackApprox Approximation to dividend option pricing
%   CODE WRITTEN BY MAITREYI MANDAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the dividend dates 
div=D{2,1}; % current dividend in $
ctr=1;
dte(1)=D{2,2};dte(2)=dte(1)+D{2,3};dte(3)=dte(2)+D{2,3};
dte=dte/12;
if dte(3)<T ctr=ctr+3; dte(4)=T;
elseif dte(3)>T && dte(2)<T ctr=ctr+2; dte(3)=T;
else ctr=ctr+1; dte(2)=T;
end

% discounted dividend
dDiv_tn=0;
if ctr >2 dDiv_tn = sum(div.*exp(-r.*(dte(1:ctr-2)))); end% for t(n) - last ex-dividend date
S0tn=S0-dDiv_tn;
dDiv_T = sum(div.*exp(-r.*(dte(1:ctr-1)))); % for T - maturity date
S0T=S0-dDiv_T;

% calculate the call Option prices on dates derives above
opttn=blsprice(S0tn,K,r,dte(ctr-1),si);
optT=blsprice(S0T,K,r,dte(ctr),si);
% end
option=max(opttn,optT);

end

