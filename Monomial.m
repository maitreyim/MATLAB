function [ M ] = Monomial(x,n)
M=nan(size(x,1),n);
if(n>0)
    M(:,1)=x.*x - x.*x + 1;
end

if(n>1)
    M(:,2)=x.^1;
end

if(n>2)
    M(:,3)=x.^2;
end
if(n>3)
    M(:,4)=x.^3;
end
if(n>4)
    M(:,5)=x.^4;
end
end

    
