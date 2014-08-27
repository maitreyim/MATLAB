function [u,d,p] = Q1c_bino(sigma,r,delta)
u=exp((r-(sigma^2)/2)*delta+sigma*sqrt(delta));
d=exp((r-(sigma^2)/2)*delta-sigma*sqrt(delta));
p=0.5;
end

