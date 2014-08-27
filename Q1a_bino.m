function [u,d,p] = Q1a_bino(sigma,r,delta)
c=0.5*(exp(-r*delta)+exp((r+sigma^2)*delta));
d=c-sqrt(c^2-1);
u=1/d;
p=((exp(r*delta)-d)/(u-d));
end

