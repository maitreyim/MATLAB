function [ rmse ] = calcVasicekNSquareRMSE( x, cmt)
%calcVasicekNSquareRMSE Calculate NSquare RMSE for Interest Rates
%   Function to calculate the Root Mean Square Error for Vasicek NSquare
% method
% cmt - data for constant maturity treasury
% alphax - 
% betax - 
% sigmax - 
% betay - 
% sigmay - 

alphax = x(1); betax = x(2); sigmax = x(3); betay = x(4); sigmay = x(5);
T1 = 0.25; T2 = 10;
alphay = 0;
%STEP 2:: Compute Ax(0.25) Ay(0.25) Bx(0.25) By(0.25)
ATf = AT(); BTf = BT();
AxT1 = ATf(alphax,betax,sigmax,T1);
BxT1 = BTf(betax,T1);
AyT1 = ATf(alphay,betay,sigmay,T1);
ByT1 = BTf(betay,T1);

AxT2 = ATf(alphax,betax,sigmax,T2);
BxT2 = BTf(betax,T2);
AyT2 = ATf(alphay,betay,sigmay,T2);
ByT2 = BTf(betay,T2);

% STEP 2 a (used in STEP 3 and 4):: calc Ax(T), Ay(T), Bx(T), By(T)
% used to calculate D(T)
for t=0.5:0.5:10
    Ax(t*2) = ATf(alphax,betax,sigmax,t);
    Bx(t*2) = BTf(betax,t);
    Ay(t*2) = ATf(alphay,betay,sigmay,t);
    By(t*2) = BTf(betay,t);
end

% STEP 3 and 4 :: Find the time series of X and Y
% define simultaneous equations
A = [ BxT1/T1 ByT1/T1;
    BxT2/T2 ByT2/T2];
for i=1:size(cmt);

    b = [ cmt(i,1)+(log(AxT1)/T1)+(log(AyT1)/T1);
        cmt(i,6)+(log(AxT2)/T2)+(log(AyT2)/T2)];
    % solve for X and Y using simultaneous equations
    v(i,:)=(A\b)';
end

% STEP 5 :: compute par rates for T = 2 3 5 7
% Compute D(T) for t = 0.5 to 10 for the 650 rows
for i=1:size(cmt)
    D(i,:) = Ax.*Ay.*exp(-Bx.*v(i,1)-By.*v(i,2));
end
% compute coupon rates for T = 2 3 5 7
% define vector of T
T = [2 3 5 7];
for i=1:size(cmt)
    for j=1:4
        columnForD = T(j)*2;
        cumsumD = cumsum(D(i,1:columnForD));
        C(i,j) = (2*(100 - 100.*D(i,T(j)))/cumsumD(columnForD))*(1/100);
    end
end

%STEP 6 :: compute RMSE for each day and for all dates
MSE = (cmt(:,2) - C(:,1)).^2 + (cmt(:,3) - C(:,2)).^2 + ...
        (cmt(:,4) - C(:,3)).^2 + (cmt(:,5) - C(:,4)).^2;
rmse = sqrt(mean(MSE));    
end

