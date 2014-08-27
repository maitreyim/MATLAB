Xbar = 20;
X0 = 10;
K = 20;
T = 1;
r = .04;
sigma = .25;
% blsprice(X0, K, r, T, sigma)
% 0.004841559000219
 
 
% M = Xbar;
% N = T/.002 - 1;

% M = 500;
% N = 500;


dx = Xbar/(1+M);
dt = T/(1+N);

xgrid = linspace(0, Xbar, M+1)';
tgrid = linspace(0, T   , N+1)';

a = -xgrid*r/(4*dx) + (xgrid * sigma).^2/(4*dx^2);
c =  xgrid*r/(4*dx) + (xgrid * sigma).^2/(4*dx^2);
b = -r/2      - (xgrid * sigma).^2/(2*dx^2) - 1/dt;
d =  r/2      + (xgrid * sigma).^2/(2*dx^2) - 1/dt;

A = diag(b) + diag(c(1:(end-1)),1) + diag(a(2:end),-1);
A(end, end-1) = 0;
A(end, end) = 1;

A(1,1) = 1;
A(1,2) = 0;

F = zeros(M+1, N+1);
F(:, end) = max(0, xgrid - K);
endval = (a(end) + b(end)) * (Xbar - K);
% F(:, end) = max(K-xgrid,0);
% endval = (a(end) + b(end)) * (K-Xbar);

for n = N:(-1):1 
  Dn = [0 ; ...
         -a(2:(end-1),1) .* F(1:(end-2),n+1) + ...
          d(2:(end-1),1) .* F(2:(end-1),n+1) + ...
         -c(2:(end-1),1) .* F(3:(end  ),n+1); ...
         Xbar - K + dx];

  F(:,n) = A\Dn;
  if( 0 == mod( n, 10 ) ) 
      disp(['iteration ' num2str(n) ' complete']) 
  end
end

mesh(tgrid', xgrid, F)
F(X0 == xgrid ,find( 0 == tgrid))


