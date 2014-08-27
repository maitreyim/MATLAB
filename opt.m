[num,txt,raw]=xlsread('C:\Users\maitreyi\Desktop\opt.xlsx');

vol=num(:,5);
rate=num(:,4);
strike=num(:,2);
price=num(:,1);
time=num(:,3);
for i=1:1261
c(i)=blsprice(price(i),strike(i),rate(i),time(i),vol(i));
end
c;

