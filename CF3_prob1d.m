% Problem 1, Part d
%Compute P(Y2>5)
T=2;%Final Time
N=1000;%# of obs
Y0=3/4;% Initial value of Y
delta=T/N;% Creating delta
path=2000;% # of paths
values=nan(path,1);
for j=1:path;% loop over path
z=randn(N,1);%creating random normal numbers
Y = zeros(N+1,1);
Y(1)=Y0;
for i=1:N % time loop
    Y(i+1,1)=Y(i,1)+((2/(1+(i-1)*delta))*Y(i,1)+ (1+(i-1)*delta^3)/3)*delta+(1+((i-1)*(delta^3))/3)*sqrt(delta)*z(i,1);
    % SDE of Yt as given in the question
end
values(j,1)=Y(N,1);
end;
E5=mean(values>5);% Computing Probability
