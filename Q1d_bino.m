function [u,d,p] = Q1d_bino(sigma,r,delta)
u=exp(sigma*sqrt(delta));
d=exp(sigma*sqrt(delta)*(-1));
p= 0.5+0.5*(((r-(sigma^2)/2)*sqrt(delta))/sigma);
end

