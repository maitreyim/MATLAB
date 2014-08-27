%calculates fair price of variance swap for give input volatilities
%Based on paper "More Than You Ever Wanted To Know About Volatility Swaps" by Derman et al
%Implements example of table 1:
%The portfolio of European-style put and call options used for
%calculating the cost of capturing realized variance in the presence 
%of the implied volatility skew with a discrete set of options strikes.
%For more details refer
%http://www.quantcode.com/modules/mydownloads/singlefile.php?cid=11&lid=72

%Note : this file needs function BlackScholesPrice which can be downloaded from
%http://www.quantcode.com/uploads/BlackScholesPrice.m

%-----------------------Enter parameters
S0=100; %spot price

putstrikesvec=100:-5:45;  %enter vector for available strikes for puts
volvec_put=0.2:0.01:0.3; %enter vector of implied vols for puts

callstrikesvec=100:5:140; %enter vector for available strikes for calls
volvec_call=.2:-0.01:0.13; %enter vector of implied vols for calls

r=0.05; %risk free rate
T=90/365; %maturity of 3 months
SQ=100; %strike price which is nearest to forward price
%--------------------------------------------

callstrikesvec=callstrikesvec';
putstrikesvec=putstrikesvec';

nc=size(callstrikesvec,1)-1;
np=size(putstrikesvec,1)-1;
%create a display matrix to store values of Table 1
resultsmat=zeros(np+nc,5);

%-----Creating portfolio of call options
%use equation (A4) on page 42 to find function value
fvec=callstrikesvec*0;

wck=fvec(1:size(fvec)-1); %weights of calls

ST=callstrikesvec(1);
fvec(1)=(2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) ;

for i=2:size(callstrikesvec)
  ST=callstrikesvec(i);
  
  fvec(i)=(2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) ;
  %use eq (A7) on page 43
  wck(i-1)=(fvec(i)-fvec(i-1))/(callstrikesvec(i)-callstrikesvec(i-1));
  if i>2
    wck(i-1)=wck(i-1)-sum(wck(1:i-2));
  end
  end

counter=np+1;

call_cost=0; %total value of call portfolio
CallPutFlag=c;
for i=1:size(wck)
  v=volvec_call(i);
  X=callstrikesvec(i);
  call_price=BlackScholesPrice(CallPutFlag,S0,X,T,r,v);
  call_cost=call_cost+call_price*wck(i);
  resultsmat(counter,1)=X;
  resultsmat(counter,2)=v;
  resultsmat(counter,3)=wck(i)*10000;
  resultsmat(counter,4)=call_price;
  resultsmat(counter,5)=call_price*wck(i)*10000;
  counter=counter+1;
end


%-----Creating portfolio of put options
%use equation (A4) on page 42 to find function value
fvec=putstrikesvec*0;

wpk=fvec(1:size(fvec)-1); %weights of calls

ST=putstrikesvec(1);
fvec(1)=(2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) ;

for i=2:size(putstrikesvec)
  ST=putstrikesvec(i);
  fvec(i)=(2/T)*( (ST-SQ)/SQ - log(ST/SQ) ) ;
  %use eq (A8) on page 43
  wpk(i-1)=(fvec(i)-fvec(i-1))/(putstrikesvec(i-1)-putstrikesvec(i));
  if i>2,
    wpk(i-1)=wpk(i-1)-sum(wpk(1:i-2));
    end
end

counter=np;

put_cost=0; %total value of put portfolio
CallPutFlag=p;
for i=1:size(wpk),
  v=volvec_put(i);
  X=putstrikesvec(i);
  put_price=BlackScholesPrice(CallPutFlag,S0,X,T,r,v);
  put_cost=put_cost+put_price*wpk(i);
  
  resultsmat(counter,1)=X;
  resultsmat(counter,2)=v;
  resultsmat(counter,3)=wpk(i)*10000;
  resultsmat(counter,4)=put_price;
  resultsmat(counter,5)=put_price*wpk(i)*10000;
  counter=counter-1;
end


%use eqn. (29) to find total weighted cost of portfolio
postfolio_cost=put_cost+call_cost;

%verify table 1 values are reproduced
resultsmat
Total_Cost=sum(resultsmat(:,5))

%use eqn (27) to find fair rate for variance swap
Kvar=(2/T)* ( r*T-(S0*exp(r*T)/SQ-1)-log(SQ/S0) ) + exp(r*T)*(postfolio_cost) 
rough_fairvol=Kvar^0.5 %rough estimate of fair volatility



