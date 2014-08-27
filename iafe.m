[NUMERIC,TXT,RAW]=xlsread('C:\Users\maitreyi\Desktop\iafe_new.xlsx','Sheet1','a1:e5');
corr=cell2mat(RAW);

EigenV= eig(corr);