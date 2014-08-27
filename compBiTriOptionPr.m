function [ cpPrice, stockVec, cvVec, evVec, evcvVec ] = compBiTriOptionPr( r, si, S0, K, T, u, d, pu, pd, pm, n, treeTyp, optTyp, CallOrPut )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if treeTyp == 'B' fStep = (n+1)*(n+1+1)/2;
else if treeTyp == 'T' fStep=(n+1)^2; 
    end
end
l=0;
h=1;
t=T/n;

for i=0:n
    if treeTyp == 'B' 
        j =i+1; end
    if treeTyp == 'T' 
        j = 2*i + 1;
    end
    for k=j:-1:1
        l=l+1;
        mid = j/2;
        ctrT = (j-k+1);
        if treeTyp == 'B'
            powu = k-1;
            powd = j-powu-1;
        end
        % look at trinomial logic later
        if treeTyp == 'T'
            if ceil(mid-ctrT) < 0
                powu = 0;
            else
                powu = ceil(mid-ctrT);
            end;
            if floor(ctrT-mid) < 0
                powd = 0;
            else
                powd = floor(ctrT-mid);
            end;
        end
        S(l) = S0*(u^powu)*(d^powd);
    end
end
m = fStep;

for i=n:-1:0
    if treeTyp == 'B' 
        j =i+1; end
    if treeTyp == 'T' 
        j = 2*i + 1; end     
    m = m-j;
    l = m;
    for k=1:j
        % calculate EV and CV
        l = l + 1;
        % Call and put Value
        if CallOrPut == 'C' 
            if S(l) > K
                evVec(l) = S(l) - K;
            else
                evVec(l) = 0;
            end
        else 
            if S(l) < K
                evVec(l) = K - S(l);
            else
                evVec(l) = 0;
            end
        end
        %call and put value ends
        %set cv = 0 if you are in final step
        if i == n
            cvVec(l) = 0;
            evcvVec(l) = evVec(l);
        else
            if treeTyp=='B'
                cvVec(l) = exp(-r*t)*(pu*evcvVec(l+j) + pd*evcvVec(l+j+1));
            else if treeTyp=='T'
                    cvVec(l) = exp(-r*t)*(pu*evcvVec(l+j) + pm*evcvVec(l+j+1) + pd*evcvVec(l+j+2));
                end
            end
            if optTyp == 'E'
                evVec(l) = 0;
            end            
            if cvVec(l) > evVec(l)
                evcvVec(l) = cvVec(l);
            else
                evcvVec(l) = evVec(l);
            end

        end
            

    end
end
cpPrice = evcvVec(1);
stockVec = S;
%cvVec = zeros(fStep);
%evVec = zeros(fStep);
end

