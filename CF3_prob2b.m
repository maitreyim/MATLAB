%Question 2b
% Compute the expected value of Yt=E(1+y3)^1/3
rng(1);
W=sqrt(3)*randn(1000,1);% Generating First Sequence of Standard Wiener Process
rng(2);
Z=sqrt(3)*randn(1000,1);
Y= exp(-0.08*3+(1/3)*W+(3/4)*Z); % Generating Second Sequence of Standard Wiener Process
S=mean(Y); % Computing the Expected Value



