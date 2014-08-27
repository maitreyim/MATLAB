subplot(1,1,1);
n = 1000;
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho, n);
U = normcdf(Z,0,1);
X = [gaminv(U(:,1),2,1) betainv(U(:,2),2,2) tinv(U(:,3),5)];
plot3(X(:,1),X(:,2),X(:,3),'.');
grid on; view([-55, 15]);
xlabel('U1'); ylabel('U2'); zlabel('U3');

Sigma = [ 0.002838946 , -0.000544013, -4.7455E-05; ...
         -0.000544013,  0.022700991,  0.000282068; ...
         -4.7455E-05 ,  0.000282068,  0.000368527 ];

     
data = xlsread('C:\Users\maitreyi\Desktop\copula\data_copula.xlsx');


		
data = csvread('C:\Users\maitreyi\Desktop\copula\data_copula.csv');
Sigma = cov(data);


weights = [1,1,1]/3;

high20 = [0.1180, 0.1240,0.0920,0.1190,0.1400,0.1160,0.1240,0.1160,0.1070,0.0890,0.1130,0.1020,0.1160,0.1200,0.1160,0.1340,0.1040,0.059, 0.031];
bond = [5.87, 7.09, 6.57, 6.44, 6.35, 5.26, 5.65, 6.03, 5.02, 4.61, 4.01, 4.27, 4.29, 4.80, 4.63, 3.66, 3.26, 3.22, 2.78]/100;
sp = [.086253758, 0.019507983, 0.176577907, 0.237724286, 0.302654735, 0.242801369, 0.222782128, 0.075256342, -0.163282465, -0.167679914, -0.028885043, 0.171378842, 0.067730951, 0.085509803, 0.127230133, -0.174080518, -0.222935314, 0.20243658, 0.111994175];

sample_covariance = cov([ high20', bond', sp' ]);
sample_mean = mean([ high20', bond', sp' ],1);

portfolio_sd = weights * sample_covariance * weights';
portfolio_mean = sample_mean * weights';
var95 = portfolio_mean - 1.65 * portfolio_sd;


scaled_high20 =( high20 - min(high20))/(max(high20) - min(high20));

alpha1 = 2.2458;
alpha2 = 0.74438;

lower = 2.5406/100;
upper = 7.3294/100;

bond_quantile = (bond - lower)/(upper - lower);
high20_quantile = betainv(scaled_high20, alpha1, alpha2);

covar_matrix = cov( [high20' bond']);
bond_normal = mvnorminv(bond, bond_mean, bond_var);




