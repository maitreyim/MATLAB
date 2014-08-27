% PROJECT 4
% main script for calling Project 4 functions
% call option prcing function
%[ cpPrice, stockVec, cvVec, evVec  ] = compBiTriOptionPr( r, si, S0, K, T, u, d, pu, pd, pm, n, treeTyp, optTyp, CallOrPut )
clear all;
T = 0.5; r = 0.05; si=0.24; S0=32; K=30;
n10 = 10; n20=20; n40=40; n80=80; n100=100; n200=200; n500=500;
q1nc = zeros(7, 5);
q1nc(:,1) = [10, 20, 40, 80, 100, 200, 500];
disp(q1nc);
treeTyp='B'; optTyp='E'; CallOrPut='C';
for i=1:7
    [u,d,p] = Q1a_bino(si,r,T/q1nc(i,1));
    [q1nc(i,2), stockVec, cvVec, evVec, evcvVec  ] = compBiTriOptionPr( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
   % error1(i)=[q1nc(i,2), stockVec, cvVec, evVec, evcvVec  ]-blsprice(32,30,T/q1nc(i,1);
    
    [u,d,p] = Q1b_bino(si,r,T/q1nc(i,1));
    [q1nc(i,3), stockVec, cvVec, ~, evcvVec  ] = compBiTriOptionPr( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
     
    [u,d,p] = Q1c_bino(si,r,T/q1nc(i,1));
    [q1nc(i,4), stockVec, cvVec, evVec, evcvVec  ] = compBiTriOptionPr( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
     
    [u,d,p] = Q1d_bino(si,r,T/q1nc(i,1));
    [q1nc(i,5), stockVec, cvVec, evVec, evcvVec  ] = compBiTriOptionPr( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
end
disp('Stock Vector is : '); 
disp(stockVec);
disp('Current Value vector is : ');
disp(cvVec);
disp('Expected Value vector is : ');
disp(evVec);

disp('EvCv vector is : ');
disp(evcvVec);

disp('Call Price is : ');
disp(q1nc);

figure(1);
plot(q1nc(:,1),q1nc(:,2),q1nc(:,1),q1nc(:,3),q1nc(:,1),q1nc(:,4),q1nc(:,1),q1nc(:,5));

