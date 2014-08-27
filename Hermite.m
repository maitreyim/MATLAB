function [ H ] = Hermite(x,n)
H=nan(size(x,1),n);
if(n>0)
    H(:,1)=x.*x - x.*x +1;
end

if (n>1)
    H(:,2)=2.*x;
end

if(n>2)
    H(:,3)=(4.*x.^2) - 2;
end
if (n>3)
    H(:,4)=8.*x.^3 - 12.*x;
end

if(n>4)
    H(:,5)=16.*x.^4 - 56.*x.^2 + 16;
end
   
 end