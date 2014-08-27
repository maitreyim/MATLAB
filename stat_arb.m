%HW1

%read presorted datastream and bblist
[num, txt, raw]=xlsread('DataStream.xlsx','Sheet1');
[num2, txt2, raw2]=xlsread('B2.xlsx','Sheet1');
[I J]=size(raw);
for i=1:I
for j=1:J
raw{i,j}=num2str(raw{i,j});
end
end
dscode=unique(raw(2:end,3));
allstocks=struct('dscode',dscode);

for i=1:length(dscode)
allstocks(i).dscode=dscode(i);
end

i2=0;

%insert date and so on
for n=2:length(raw)
i=strmatch(raw{n,3},dscode);
if i~=i2;
j=1;
allstocks(i).namelist(1).date=datestr(raw{n,1},'dd-mmm-yyyy');
allstocks(i).namelist(1).name=raw{n,4};
allstocks(i).industrylist=raw{n,7};
allstocks(i).ibeslist=raw{n,5};
allstocks(i).isinlist=raw{n,6};
allstocks(i).indexlist(1).index=raw(n,10);
elseif i==i2 & strcmp(raw{n,4},raw{n-1,4})==0;
j=j+1;
allstocks(i).namelist(j).name=raw{n,4};
allstocks(i).namelist(j).date=datestr(raw{n,1},'dd-mmm-yyyy');
end
end

%insert bblist
for n=2:length(raw2)
if strcmp(num2str(raw2{n,1}),'')~=1
i=strmatch(num2str(raw2{n,1}),dscode);
allstocks(i).bblist=num2str(raw2{n,2});
end
end