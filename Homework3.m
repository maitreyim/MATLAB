%Computational Finance MGMT 237G
%Homework 3

%Question 3, Part c)
call = blsprice(90.0, 105.0, 0.04, 5.0, 0.2);   
epsilon = 0.0005;
r = 0.04;
sigma = 0.25;
X = 20;
T = 0.5;
delta = zeros(1, 11);
gamma = zeros(1, 11);
theta = zeros(1, 11);
vega = zeros(1, 11);
rho = zeros(1, 11);

deltaCom = zeros(1, 11);
gammaCom = zeros(1, 11);
thetaCom = zeros(1, 11);
vegaCom = zeros(1, 11);
rhoCom = zeros(1, 11);

for i = 15:1:25
    S0 = i;

    %Delta
    S1 = S0;
    S2 = S0 + epsilon;
    CofS1 = blsprice(S1, X, r, T, sigma);
    CofS2 = blsprice(S2, X, r, T, sigma); 
    delta1 = (CofS2 - CofS1) / epsilon;
    delta(i - 14) = delta1;
    
    d1 = (log(S0/X) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T));
	d2 = d1 - sigma * sqrt(T);
    
    N1 = ApproximationOfNormalDistribution(d1);
    N2 = ApproximationOfNormalDistribution(d2);
    deltaCom(i - 14) = N1;
    %Gamma
    S3 = S0 + 2 * epsilon;
    CofS3 = blsprice(S3, X, r, T, sigma);
    delta2 = (CofS3 - CofS2) / epsilon;
    gamma(i - 14) = (delta2 - delta1) / epsilon;
    
    gammaCom(i - 14) = N1 / (S0*sigma*sqrt(T));
    %Theta
    T1 = T;
    T2 = T + epsilon;
    CofT1 = blsprice(S0, X, r, T1, sigma);
    CofT2 = blsprice(S0, X, r, T2, sigma); 
    theta(i - 14) = (CofT2 - CofT1) / epsilon;

    thetaCom(i - 14) = -S0*sigma*N1/(2*sqrt(T)) - r*X*exp(-r*T)*N2;
    %Vega
    sigma1 = sigma;
    sigma2 = sigma + epsilon;
    CofSigma1 = blsprice(S0, X, r, T, sigma1);
    CofSigma2 = blsprice(S0, X, r, T, sigma2); 
    vega(i - 14) = (CofSigma2 - CofSigma1) / epsilon;
    
    vegaCom(i - 14) = S0 * sqrt(T)*N1;
    %Rho
    r1 = r;
    r2 = r + epsilon;
    CofR1 = blsprice(S0, X, r1, T, sigma);
    CofR2 = blsprice(S0, X, r2, T, sigma); 
    rho(i - 14) = (CofR2 - CofR1) / epsilon;
    
    rhoCom(i - 14) = X*T*exp(-r*T)*N2;
end;

figure(1)
plot(15:1:25, delta, '-r', 15:1:25, gamma, '-b', 15:1:25, theta, '-g', 15:1:25, vega, '-k', 15:1:25, rho, '-c', 'linewidth',2);
legend('Delta', 'Gamma', 'Theta', 'Vega', 'Rho');
figure(5)
plot(15:1:25, deltaCom, '-r', 15:1:25, gammaCom, '-b', 15:1:25, thetaCom, '-g', 15:1:25, vegaCom, '-k', 15:1:25, rhoCom, '-c', 'linewidth',2);
legend('Delta', 'Gamma', 'Theta', 'Vega', 'Rho');

figure(6)
plot(15:1:25, gamma, '-r', 15:1:25, gammaCom, '-b', 'linewidth',2);

%Question 5
numberOfVectors = 100;
%Part a)
vectorOfUniforms = zeros(numberOfVectors, 2);
seed = 98242;
iterations = 200;
uniforms = zeros(1, iterations);
uniforms = LGM(iterations, seed);
count = 1;

for i = 1:numberOfVectors
    vectorOfUniforms(i, 1) = uniforms(count);
    count = count + 1;
    vectorOfUniforms(i, 2) = uniforms(count);
    count = count + 1;
end;

%Part b)
haltonSeqB = zeros(numberOfVectors, 2);
seqBase2 = GetHalton(numberOfVectors, 2);
seqBase7 = GetHalton(numberOfVectors, 7);

for i = 1:numberOfVectors
    haltonSeqB(i, 1) = seqBase2(i);
    haltonSeqB(i, 2) = seqBase7(i);
end;

%Part c)
%Part b)
haltonSeqC = zeros(numberOfVectors, 2);
seqBase2 = GetHalton(numberOfVectors, 2);
seqBase4 = GetHalton(numberOfVectors, 4);

for i = 1:numberOfVectors
    haltonSeqC(i, 1) = seqBase2(i);
    haltonSeqC(i, 2) = seqBase4(i);
end;

%figure(2)
%plot(vectorOfUniforms(:,1), vectorOfUniforms(:,2));
%figure(3)
%plot(haltonSeqB(:,1), haltonSeqB(:,2));
%figure(4)
%plot(haltonSeqC(:,1), haltonSeqC(:,2));

%Part e)
n = 10000;
baseX = 5;
baseY = 7;
seqBaseX = GetHalton(n, baseX);
seqBaseY = GetHalton(n, baseY);
fxyi = 0;
for i = 1:n
    x = seqBaseX(i);
    y = seqBaseY(i);
    sinPart = sin(6*pi*x);
    cosPart = nthroot(cos(2*pi*y),3);
    fxyi = fxyi + exp(-x*y) * (sinPart + cosPart);
end;

integral = fxyi / n;
dd = 0;

    
