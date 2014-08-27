%CDS pricing function
%a=Beginning time of the CDS (usually 0)
%b=Maturity in years of the CDS
%s=Market quote
%L=Loss given default
%gamma=Intensity vector
%PEURC=Default-free zero-coupon bond prices
function z=CDS(a,b,s,L,gamma,PEURC);
j=(b-a)*4;
Gamma_j=zeros(1,j+1);
for i=2:(j+1)
Gamma_j(i)=Gamma_j(i-1)+gamma(i);
end
fir=0;
sec=0;
thir=0;
for i=1:j
fir=fir+gamma(i)*quad(@(u)FIRST(u,Gamma_j,gamma,PEURC,i-1),i-1,i);
sec=sec+exp(-Gamma_j(i))*PEURC(i*90);
thir=thir+gamma(i)*quad(@(u)THIRD(u,Gamma_j,gamma,PEURC,i-1),i-1,i);
end
fir=s*fir;
sec=s*sec;
thir=L*thir;
z=fir+sec-thir;
return;
%First argument auxiliary function for the CDS pricing
function y=FIRST(u,Gamma_j,gamma,PEURC,mmin)
days=floor(u*90);
days(days == 0) = 1;
y=exp(-Gamma_j(mmin+1)-gamma(mmin+2).*(u-mmin)).*PEURC(days)'.*(u-mmin);
return;
%Third argument auxiliary function for the CDS pricing
function y=THIRD(u,Gamma_j,gamma,PEURC,mmin)
days=floor(u*90);
days(days == 0) = 1;
y=exp(-Gamma_j(mmin+1)-gamma(mmin+2).*(u-mmin)).*PEURC(days)';
return;
%Code for stripping out constant intensities
%Function to strip out constant intensities
%L=Loss given default
%ba=Maturity terms
%sa=Market quotes
%gammai=Intensity parameters
function z=calib(PEURC);
L=1-0.15;
ba=[1 3 5 7 10];
%sa=[.01925 .0215 .0225 .0235 .0235];
%sa=[.0725 .0630 .0570 .0570 .0570];
%sa=[.1450 .1200 .0940 .0850 .0850];
sa=[.5050 .2100 .1500 .1250 .1100];
gamma1=.05;
gamma2=.05;
gamma3=.05;
gamma4=.05;
gamma5=.05;
gamma1=fzero(@(gamma1)CDS(0,ba(1),sa(1)/4,L,[[gamma1/4*ones(1,5)]
[gamma2/4*ones(1,8)] [gamma3/4*ones(1,8)] [gamma4/4*ones(1,8)]
[gamma5/4*ones(1,12)]],PEURC),gamma1);
gamma2=fzero(@(gamma2)CDS(0,ba(2),sa(2)/4,L,[[gamma1/4*ones(1,5)]
[gamma2/4*ones(1,8)] [gamma3/4*ones(1,8)] [gamma4/4*ones(1,8)]
[gamma5/4*ones(1,12)]],PEURC),gamma2);
gamma3=fzero(@(gamma3)CDS(0,ba(3),sa(3)/4,L,[[gamma1/4*ones(1,5)]
[gamma2/4*ones(1,8)] [gamma3/4*ones(1,8)] [gamma4/4*ones(1,8)]
[gamma5/4*ones(1,12)]],PEURC),gamma3);
gamma4=fzero(@(gamma4)CDS(0,ba(4),sa(4)/4,L,[[gamma1/4*ones(1,5)]
[gamma2/4*ones(1,8)] [gamma3/4*ones(1,8)] [gamma4/4*ones(1,8)]
[gamma5/4*ones(1,12)]],PEURC),gamma4);
gamma5=fzero(@(gamma5)CDS(0,ba(5),sa(5)/4,L,[[gamma1/4*ones(1,5)]
[gamma2/4*ones(1,8)] [gamma3/4*ones(1,8)] [gamma4/4*ones(1,8)]
[gamma5/4*ones(1,12)]],PEURC),gamma5);
z=[gamma1;gamma2;gamma3;gamma4;gamma5];
return;