[NUMERIC,TXT,RAW]=xlsread('C:\Users\maitreyi\Desktop\FX.xls','Ret_D','b2:e732');
Assetret=cell2mat(RAW);
y1=Assetret(:,1);
y2=Assetret(:,2);
y3=Assetret(:,3);
y4=Assetret(:,4);
T1=length(y1);
% MA Estimation
W=20;
for t=T1-5:T1
	t1=t-W+1;
	window=y1(t1:t);	%  estimation window
	
	
end
sigma=std(window);
% EWMA Estimation
lambda = 0.94;	
s11 = var(y1(1:30));	% initial variance
for t = 2:T1							
	
end
s11 = lambda * s11  + (1-lambda) * y1(t-1)^2;
