function [Xu,Xd,pu,pd,pm] = Q6b_trino(sigma,r,delta)

Xu=sigma*sqrt(3*delta);
Xd= exp(-sigma*sqrt(3*delta));
pd=0.5*(((sigma^2)*delta+((r-(sigma^2)/2)^2)/(Xu^2))-((r-(sigma^2)/2)*delta)/Xd);
pu=0.5*(((sigma^2)*delta+((r-(sigma^2)/2)^2)/(Xu^2))+((r-(sigma^2)/2)*delta)/Xd);
pm=(1-pu-pd);

end

