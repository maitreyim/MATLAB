%Question 5
%Part A: Generation of 2-dimensional Uniform Numbers
N_vec = 100; % Length of vector
%Part a)
v_Uniforms = zeros(N_vec , 2); % Generate Uniforms
seed = 98242;% Use the seed
trials = 200;% # of iterations
uniforms = zeros(1, trials);% Generate Uniforms
uniforms = LGM(iterations, seed);% Finally Generate Uniforms
count = 1;
% For loop 
for i = 1:N_vec 
    v_Uniforms(i, 1) = uniforms(count);
    count = count + 1;
    v_Uniforms(i, 2) = uniforms(count);
    count = count + 1;
end;
%Part b: Generation of Halton Sequence of base 2 & 7
haltonSeqB = zeros(N_vec , 2);% Create The Vector
seqBase2 = GetHalton(N_vec , 2); % Use the function GetHalton (base 2) as defined
seqBase7 = GetHalton(N_vec , 7);% Use the function GetHalton (base 7) as defined

for i = 1:N_vec
    haltonSeqB(i, 1) = seqBase2(i);
    haltonSeqB(i, 2) = seqBase7(i);
end;
%Part c: Generation of Halton Sequence of base 2 & 4
haltonSeqC = zeros(N_vec , 2);% Create The Vector
seqBase2 = GetHalton(N_vec , 2);% Use the function GetHalton (base 2) as defined
seqBase4 = GetHalton(N_vec , 4);% Use the function GetHalton (base 7) as defined

for i = 1:N_vec
    haltonSeqC(i, 1) = seqBase2(i);
    haltonSeqC(i, 2) = seqBase4(i);
end;

% Plotting them all
subplot(2,2,1) % Using Subplot to plot them all in one graph
scatter(vectorOfUniforms(:,1), vectorOfUniforms(:,2));% Creating Scatter Diagram
Title('Scatterplot of Uniform Numbers');% Title

subplot(2,2,2) % Using Subplot to plot them all in one graph
scatter(haltonSeqB(:,1), haltonSeqB(:,2)) % Creating Scatter Diagram
Title('Halton Sequence of Base 2,7');% Title

subplot(2,2,3)% Using Subplot to plot them all in one graph
scatter(haltonSeqC(:,1), haltonSeqC(:,2))% Creating Scatter Diagram
Title('Halton Sequence of Base 2,4');% Title

%Part e: Compute the Integral as asked in the question 5e with different
%touples
n = 10000;% # of obs
Touple_X = 5; % First sequence of touple (X,Y)
Touple_Y = 7; % Second  sequence of touple (X,Y)
seq_toupleX = GetHalton(n, Touple_X); % Generating Sequence for X
seq_toupleY = GetHalton(n, Touple_Y); % Generating Sequence for Y
fxyi = 0;
for i = 1:n
    x = seq_toupleX(i);
    y = seq_toupleY(i);
    sin= sin(6*pi*x);
    cos_p= nthroot(cos(2*pi*y),3);
    fxyi = fxyi + exp(-x*y) * (sin + cos_p);
end;

Value_integral = fxyi / n;% Computing the integral
dd = 0;
