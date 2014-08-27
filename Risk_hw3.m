vixd=xlsread('C:\Users\maitreyi\Desktop\FX.xls','Ret_D','b1:e733');
vixw=xlsread('C:\Users\maitreyi\Desktop\FX.xls','Ret_W','b1:e260');
vixm=xlsread('C:\Users\maitreyi\Desktop\FX.xls','Ret_M','b1:e59');

spec = garchset('VarianceModel', 'GARCH', 'P', 1, 'Q', 1, 'MaxFunEvals', 2000, 'MaxIter', 2000);
[Coeff, Errors, LLF, Inno, Sigmas, Summ] = garchfit(spec, vixd(:,1));