% Question 1: Compare different Binomials using different u,d,p parameters
% main script for calling functions
% call option prcing function
clear all;
T = 0.5; r = 0.05; si=0.24; S0=32; K=30;
q1nc = zeros(7, 5);
q1nc(:,1) = [10, 20, 40, 80, 100, 200, 500];
%disp(q1nc);
treeTyp='B'; optTyp='E'; CallOrPut='C';
True_P =blsprice(32,30,0.05,0.5,0.24);
for i= 1:7
    
    [u,d,p] = Q1a_bino(si,r,T/q1nc(i,1));
    [q1nc(i,2), stockVec, cvVec, evVec, evcvVec  ] = Cf4_q1bino( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
    %error1(i)= repmat(True_P, size([q1nc(i,2), stockVec, cvVec, evVec, evcvVec  ]));

    
    [u,d,p] = Q1b_bino(si,r,T/q1nc(i,1));
    [q1nc(i,3), stockVec, cvVec, evVec, evcvVec  ] = Cf4_q1bino( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
        %error2(i)= repmat(True_P, size([q1nc(i,3), stockVec, cvVec, evVec, evcvVec  ]));
     
    [u,d,p] = Q1c_bino(si,r,T/q1nc(i,1));
    [q1nc(i,4), stockVec, cvVec, evVec, evcvVec  ] = Cf4_q1bino( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
     %error3(i)= repmat(True_P, size([q1nc(i,4), stockVec, cvVec, evVec, evcvVec  ]));
     
    [u,d,p] = Q1d_bino(si,r,T/q1nc(i,1));
    [q1nc(i,5), stockVec, cvVec, evVec, evcvVec  ] = Cf4_q1bino( r, si, S0, K, T, u, d, p, 1-p, 0, q1nc(i,1), treeTyp, optTyp, CallOrPut);
    %error4(i)= repmat(True_P, size([q1nc(i,5), stockVec, cvVec, evVec, evcvVec  ]));
end
% Plotting the Figure to compare different Binomial Values with u,d,p
% parameters
figure(1);
plot(q1nc(:,1),q1nc(:,2),q1nc(:,1),q1nc(:,3),q1nc(:,1),q1nc(:,4),q1nc(:,1),q1nc(:,5));
Title('Comparison of Binomials with different u,d,p Parameters');
Xlabel('# of Steps');
Ylabel('Value of European Call Using Binomial');
