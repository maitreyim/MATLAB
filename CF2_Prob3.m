%Problem 3: ESTIMATION OF THE FOLLOWING VALUES BY MONTE CARLO SIMULATION
%Part A: Find E(W(5)^2+Sin(W5))
% Standard Wiener Process Generation
rng(1);
N=5000000;
z=randn(N,1);           % set the state of randn
%Creation of Antithetic Variable z1 to use in Antithetic Variates Method
z1=(-1)*z;
m=sqrt(5)*z;
m1=sqrt(5)*z1;
a=mean(m.^2+sin(m));
%Antithetic Variates Method
a1=mean(m1.^2+sin(m1));
av1=(a+a1)/2;
%Using Control Variates Method
cv = m.^2;
c4=corr(v1, cv);
av2=mean(cv);
% Part B: Find E(e^t/2*cos(W(t))
%T=0.5
T1 = 0.5; 
W1=sqrt(T1)*z;
b1=mean((exp(T1/2)*cos(W1)));
%Antithetic Variates Variance Reduction Technique with T=0.5
W11=sqrt(T1)*z*(-1);
b11=mean((exp(T1/2)*cos(W11)));
bav1=(b1+b11)/2;
%Control Variates Method
v2 =  exp(T1/2)*W1.^2;
v1 =  exp(T1/2)*cos(W1);
c1=corr(v1, v2);
s1=(mean(v2)+b1)/2;
%T=4
T2 = 4; 
W2=sqrt(T2)*z;
b2=mean((exp(T2/2)*cos(W2)));
%Antithetic Variates Variance Reduction Technique with T=4
W12=sqrt(T2)*z1;
b12=mean((exp(T2/2)*cos(W12)));
bav2=(b2+b12)/2;
%Control Variates Method
v22 =  exp(T2/2)*W2.^2;
v12 =  exp(T2/2)*cos(W2);
c2=corr(v12, v22);
s2=(mean(v22)+b2)/2;
%T=8
T3 = 8; 
W3=sqrt(T3)*z;
b3=mean((exp(T3/2)*cos(W3)));
%Antithetic Variates Variance Reduction Technique with T=8
W13=sqrt(T3)*z1;
b13=mean((exp(T3/2)*cos(W13)));
bav3=(b3+b13)/2;
%Control Variates Technique
v33 =  exp(T3/2)*cos(W3);
v31=  exp(T3/2)*W3.^2;
c3=corr(v33, v31); 
s3=(mean(v31)+b3)/2;






