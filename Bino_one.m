function [ bino_u ] = Bino_one(u,d,c,sigma,r,p,delta)
u=1/d;
d=c-sqrt(c^2-1);
c=0.5*(exp(-r*delta))+exp(r+(sigma^2)*delta);
p=((exp(r*delta)-d)/(u-d));
end

