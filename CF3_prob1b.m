%Problem 1, Part B
%Compute the expected value of E(Y3)
T=3;%Final Time Value
N=1000;% # of obs
Y0=3/4;% Initial value of Y0 as given
delta=T/N;%creating delta
path=1000;%# of paths created
values=nan(path,1);%to contain values
for j=1:path;% loop over path
z=randn(N,1);%creating normal random variables
Y = zeros(N,1);
Y(1)=Y0;%initial value of Y as given
for i=1:N % time loop
    Y(i+1,1)=Y(i,1)+((2/(1+(i-1)*delta))*Y(i,1)+ (1+(i-1)*delta^3)/3)*delta+(1+((i-1)*(delta^3))/3)*sqrt(delta)*z(i,1);
end
    values(j,1)=Y(N,1);
end;
E2=mean(values);% Finally computing the expected value
