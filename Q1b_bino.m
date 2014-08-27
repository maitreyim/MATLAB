function [u,d,p] = Q1b_bino(sigma,r,delta)
u=exp(r*delta)*(1+sqrt(exp((sigma^2)*delta)-1));
d=exp(r*delta)*(1-sqrt(exp((sigma^2)*delta)-1));
p=0.5;
end

