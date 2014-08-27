% load('C:\Users\User-PC\Desktop\MFE\Statistical Arbitrage\HW2\allstocks.mat');
% load('C:\Users\User-PC\Desktop\MFE\Statistical Arbitrage\HW2\TickSummary.mat');
% 
%  [price, txt, raw]=xlsread('data.xlsx','UP','B6:UU1570');
% [tri, txt, raw]=xlsread('data.xlsx','RI','B6:UU1570');
% [volume, txt, raw]=xlsread('data.xlsx','UVO','B6:UU1570');
% [cap, txt, raw]=xlsread('data.xlsx','MCE','B6:UU1570');
% [MTBV, txt, raw]=xlsread('data.xlsx','MTBV','B6:UU1570');
% book=cap./MTBV;
% 
% [num, myday, raw]=xlsread('data.xlsx','myDay','B1:B1565');
% 
% [updn, txt, raw]=xlsread('Ibes codes2.xlsx','Ibes codes2','H2:H67201');
% [num, ticker, raw]=xlsread('Ibes codes2.xlsx','Ibes codes2','J2:J67201');
% [num, recdate, raw]=xlsread('Ibes codes2.xlsx','Ibes codes2','K2:K67201');
% 
% rec=zeros(1565,566);
% 
% for i=1:67200
%     for j=1:566
%         if strcmp(allstocks(j).ibeslist.ibes,ticker{i})==1
%             for k=1:1565
%                 if strcmp(myday{k},recdate{i})==1
%                     rec(k,j)=rec(k,j)+updn(i);
%                 end
%             end
%         end
%     end    
% end
% 
% spread2=mean(spread,3);
% spread3=NaN(99,566);
% for i=1:566
%     for j=1:913
%         if strcmp(dscode{j},allstocks(i).dscode)
%             spread3(:,i)=spread2(:,j);
%             break
%         end
%     end
% end
% 
% monthlist1=zeros(34,3);
% monthlist1(1:12,1)=1997;
% monthlist1(13:24,1)=1998;
% monthlist1(25:34,1)=1999;
% monthlist1(1:12,2)=(1:12)';
% monthlist1(13:24,2)=(1:12)';
% monthlist1(25:34,2)=(1:10)';
% monthlistc=[monthlist1;monthlist];
% monthlistc=monthlistc(1:72,:);
% spread4=zeros(34,566);
% for i=1:34
%     spread4(i,:)=spread3(1,:);
% end
% spread5=[spread4;spread3];
% spreadc=spread5(1:72,:);
% 
% 
% criteria_a=zeros(1565,566);
% for i=253:1565
%     criteria_a(i,:)=and((sum(isnan(price(i-252:i-1,:)))<.1*252),(sum(isnan(tri(i-252:i-1,:)))<.1*252));
% end
% criteria_b=zeros(1565,566);
% for i=253:1565
%     q=zeros(1,566);
%     for j=1:9
%         q=q+(price(i,:)==price(i-j));
%     end
%     criteria_b(i,:)=(q~=9);
% end
% criteria_c=zeros(1565,566);
% criteria_d=zeros(1565,566);
% for i=253:1565
%     criteria_c(i,:)=and((sum(isnan(volume(i-21:i-1,:)))<.1*21),(sum(isnan(tri(i-21:i-1,:)))<.1*21));
% end
% criteria_d=(~isnan(book))&(~isnan(cap));
% criteria_e=zeros(1565,566);
% for i=252:1565
%     y=year(myday{i});
%     m=month(myday{i});
%     for j=1:72
%         if y==monthlistc(j,1) && m==monthlistc(j,2)
%             criteria_e(i,:)=~isnan(spreadc(j-1,:));
%             break
%         end
%     end
% end
% criteria_f=zeros(1565,566);
% rec2=abs(rec);
% for i=253:1565
%     criteria_f(i,:)=(sum(rec2(1:i,:))>0);
% end
% 

% isactivenow1=criteria_a & criteria_b & criteria_c & criteria_e & criteria_f;

% isactivenow2=zeros(1565,566);

% monthlyactivedays=zeros(72,566);
% for i=1:1565
%     y=year(myday{i});
%     m=month(myday{i});
%     for j=1:72
%         if y==monthlistc(j,1) && m==monthlistc(j,2)
%             monthlyactivedays(j,:)=monthlyactivedays(j,:)+isactivenow1(i,:);
%             break
%         end
%     end
% end


% criteria_6a=(monthlyactivedays>0 & spreadc<0.01) | (monthlyactivedays==0 & spreadc<0.008);
% 
% count=zeros(72,566);
% dailytotal=zeros(72,566);
% m1=month(myday{1});
% j=1;
% volume2=volume;
% volume2(isnan(volume2))=0;
% for i=1:1565
%     m2=month(myday{i});
%     if m2~=m1
%         j=j+1;
%     end
%     count(j,:)=count(j,:)+(~isnan(volume(i,:)));
%     dailytotal(j,:)=dailytotal(j,:)+volume2(i,:);
%     m1=m2;
% end
% avevolume=dailytotal./count;


for i=252:1565
    y=year(myday{i});
    m=month(myday{i});
    for j=1:72
        if y==monthlistc(j,1) && m==monthlistc(j,2)
            isactivenow2(i,:)=criteria_6a(j,:) & criteria_6b(j,:);
            break
        end
    end
end

isactivenow=isactivenow1;

