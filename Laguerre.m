function [ L ] = Laguerre(x, n)

L= nan(size(x,1), n);

if( n > 0 )
  L(:, 1) = exp(-x/2);
end

if( n > 1 )
  L(:, 2) = (exp(-x/2)).*(1-x);
end

if( n > 2 )
  L(:, 3) = (exp(-x/2)).*(1-2.*x+(x.^2)/2);
end

if( n > 3 )
  L(:,4) = (exp(-x/2)).*(1-3.*x+3.*(x.^2)/2 - x.*x.*x/6);
end

if( n > 4 )
  L(:,5) = (exp(-x/2)).*(1-4.*x+3.*(x.^2) - 2.*x.*x.*x/3 + x.*x.*x.*x/24);
end

end

