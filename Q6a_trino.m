function [u,d,pu,pd,pm] = Q6a_trino(sigma,r,delta)
d=exp((-1)*sigma*sqrt(3*delta));
u=1/d;
pd=(r*delta*(1-u)+(r*delta)^2 + (sigma^2)*delta)/((u-d)*(1-d));
pu=(r*delta*(1-d)+(r*delta)^2 + (sigma^2)*delta)/((u-d)*(u-1));
pm=(1-pu-pd);

end

