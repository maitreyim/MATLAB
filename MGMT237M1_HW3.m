function MGMT237M1_HW3
% =====================================================================
% This code is created by Group 5: Jie (Edward) Sheng, Anshul Maheshwari
% Chao Liang, Li Zhang for MGMT237M1 Problem Set 3
% =====================================================================

% initialize working environment
clear all; close all; clc;
% add new working path if necessary
path(path, '')

% Step 1. Loading the Data
disp('Step 1. Initialization')
% load standard solution from PS2
load('database_v04.mat');
[T,n] = size(price);

% get index of first day of month from myday
[~,ymIdx] = unique(myday(:,4:end),'rows','first');
ymIdx = sort(ymIdx);

% find active stock index at the begin of each month
actStIdx = cell(size(ymIdx));
for i = 1:size(ymIdx,1)
    [~,actStIdx{i}] = find(isactivenow(ymIdx(i),:));
end

% calculate arithmetic return from tri
ret = NaN(T,n);
ret(2:end,:) = (tri(2:end,:)-tri(1:end-1,:))./tri(1:end-1,:);

% get industry list (P11 of Lecture 8)
ind = cell(1,n);
for i = 1:n
    ind{i} = allstocks(i).industrylist(1).industry;
end
ind = unique(char(ind),'rows');

% create industry dummy
rho = size(ind,1);
R = zeros(n,rho);
for i = 1:n
    [~,idx] = ismember(allstocks(i).industrylist(1).industry,ind,'rows');
    R(i,idx) = 1;
end

% get country list and create industry dummy
con = cell(1,n);
for i = 1:n
    con{i} = allstocks(i).indexlist(1).index;
end
con = unique(char(con),'rows');
F = zeros(n,size(con,1));
for i = 1:n
    [~,idx] = ismember(allstocks(i).indexlist(1).index,con,'rows');
    F(i,idx) = 1;
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 2. Risk Management
disp('Step 2. Risk Management')
% use myRisk() function to do calculation in this part
[sig,shrink] = myRisk(ret,actStIdx,ymIdx);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 3. Alphas
disp('Step 3. Alphas')
% use four functions to do this part
% a). short-term contrarian
alpharev = myCalAlphaRev(ret,R);

% b). short-term procyclical
alpharec = myCalAlphaRec(rec);

% c). long-term contrarian
alphaval = myCalAlphaVal(mtbv);

% d). long-term procyclical
alphamom = myCalAlphaMom(tri,ymIdx);

alphablend = 0.5.*alpharev+0.25.*alpharec+0.15.*alphaval+0.1.*alphamom;

% use myStandardize() function to d/s/w alpha
alphablend = myStandardize(alphablend);

% save work in Step 1-3, which is the Part I of this work
save('MGMT237M1_HW3_Group5_1','ymIdx','ret','R','F','actStIdx','sig',...
    'shrink','alpharev','alpharec','alphaval','alphamom','alphablend');
% load('MGMT237M1_HW3_Group5_1.mat');

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 4. Optimizer
disp('Step 4. Optimization')
% load work from Part I
% load('MGMT237M1_HW3_Group5_1');

[trade,t0,lamda,mu,booksize_op,tradesize_op] = myOptimizer(ymIdx,actStIdx,...
    sig,alphablend,volume,tcost,R,F);

% get book and trade average for evaluation of lamda and mu and long/short
bookave_op = zeros(1,3);
tradeave_op = zeros(1,3);
for i = 1:3
    idxb = (booksize_op(:,i) ~= 0);
    idxt = (tradesize_op(:,i) ~= 0);
    bookave_op(i) = mean(booksize_op(idxb,i));
    tradeave_op(i) = mean(tradesize_op(idxt,i));
end
bookave_op
tradeave_op

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 5. Backtest Results
disp('Step 5. Backtest')
[back_weight,pnl,booksize,tradesize] = myBacktest(trade,ret,tcost,myday,t0);

bookave_back = mean(booksize(booksize ~= 0))
tradeave_back = mean(tradesize(tradesize ~= 0))

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 6. Performance Statistics
disp('Step 6. Performance Statistics')
[sharpe,longest_dd,deepest_dd] = myStat(pnl,booksize,t0)

save('MGMT237M1_HW3_Group5_2','lamda','mu','t0','trade','booksize_op',...
    'tradesize_op','bookave_op','tradeave_op','back_weight','pnl',...
    'booksize','tradesize','bookave_back','tradeave_back','sharpe',...
    'longest_dd','deepest_dd');
save('MGMT237M1_HW3_Group5_submit',...
    'shrink','alpharev','alpharec','alphaval','alphamom','alphablend',...
    'lamda','mu','t0','trade','back_weight','pnl','booksize','tradesize',...
    'sharpe','longest_dd','deepest_dd');
% load('MGMT237M1_HW3_Group5_2.mat');

end

%% ********************************************************************
function [sig,shrink]= myRisk(ret,actStIdx,ymIdx)
% =====================================================================
% This function is used to calculate shrinkage estimator
% =====================================================================

sig = cell(size(ymIdx,1)-12,1);     % shrinkage estimator
shrink = zeros(size(ymIdx,1)-12,1); % shrinkage slope (beta)
for i = 13:size(ymIdx,1)    % start for the 13th month to have 1-yr past data
    if isempty(actStIdx{i}) 
        continue; 
    end
    nact = size(actStIdx{i},2);
    I = eye(nact);
    
    % get all return data of past one year for active stock 
    t = ymIdx(i)-ymIdx(i-12);
    X = ret(ymIdx(i-12):ymIdx(i)-1,actStIdx{i});
    X(isnan(X)) = 0;
    
    % get shrinkage target varbar*I (P16 of Lecture 4)
    vaI = mean(var(X)).*I;
    
    % get w^2 (P21 of Lecture 4)
    S = X'*X./t;
    XX = zeros(nact,nact,t);
    w2 = 0;
    for j = 1:t
        XX(:,:,j) = X(j,:)'*X(j,:);
        w2 = w2+norm(XX(:,:,j)-S)^2;       
    end
    w2 = w2/(t*(t-1));
    
    %estimate dispersion d^2 (P25-26 of Lecture 4)
    d2 = norm(S-vaI)^2-w2;
    %get shinkage slope (beta)
    shrink(i-12) = d2/(w2+d2);
    %get the shrinkage estimator
    sig{i-12} = (1-shrink(i-12)).*vaI+shrink(i-12).*S;
end

end

%% ********************************************************************
function alphaRev = myCalAlphaRev(ret,R)
% =====================================================================
% This function is used to calculate short-term contrarian (technical)
% d/s/w alpha
% =====================================================================

[T,n] = size(ret);
at = NaN(T,n);

% w is triangular decay weight (P13 of Lecture 8)
w = 1/11-1/231.*[0:21]';
w(end) = [];
m = eye(n)-R*inv(R'*R)*R';
% iterate for each day
for i = 21:T
    % (P15 of Lecture 8)
    ati = -sum(bsxfun(@times,w,ret(i:-1:i-20,:)))./21;
    idxnan = isnan(ati);
    at(i,~idxnan) = ati(~idxnan)*m(~idxnan,~idxnan);
end

% use myStandardize() function to d/s/w alpha
alphaRev = myStandardize(at);
end

%% ********************************************************************
function alphaRec = myCalAlphaRec(rec)
% =====================================================================
% This function is used to calculate short-term procyclical (fundamental)
% d/s/w alpha
% =====================================================================

[T,n] = size(rec);
at = zeros(T,n);
for i = 45:T
    % (P26 of Lecture 8)
    at(i,:) = sum(rec(i-44:i,:))./45;
end

% use myStandardize() function to d/s/w alpha
alphaRec = myStandardize(at);
end

%% ********************************************************************
function alphaVal = myCalAlphaVal(mtbv)
% =====================================================================
% This function is used to calculate long-term contrarian (fundamental)
% d/s/w alpha
% =====================================================================

% P7 of Lecture 5
at = 1./mtbv;
at(isinf(at)) = NaN;
% use myStandardize() function to d/s/w alpha
alphaVal = myStandardize(at);
end

%% ********************************************************************
function alphaMom = myCalAlphaMom(tri,ymIdx)
% =====================================================================
% This function is used to calculate long-term procyclical (technical)
% d/s/w alpha
% =====================================================================

[T,n] = size(tri);
at = zeros(T,n);
for i = ymIdx(13):T   % start for the 13th month to have 1-yr past data
    % get index of past 12-1 month
    m = i-ymIdx;
    [idxd,idxm] = min(m(m >= 0));
    % get return during this period
    at(i,:) = (tri(ymIdx(idxm-1)+idxd,:)-tri(ymIdx(idxm-12)+idxd,:))./...
        tri(ymIdx(idxm-12)+idxd,:);
end
at(isinf(at)) = NaN;
% use myStandardize() function to d/s/w alpha
alphaMom = myStandardize(at);
end

%% ********************************************************************
function stdAlpha = myStandardize(rawAlpha)
% =====================================================================
% This function is used to demean, standardize, windsorize alpha as showed
% on P3-4 of Lecture 5
% =====================================================================

[T,n] = size(rawAlpha);
% reiterate for better accuracy
for i = 1:5
    for j = 1:T
        % only take available alpha's, because NaN alpha can exist anywhere
        % this step can only be performed at each time
        idxnan = isnan(rawAlpha(j,:));
        if all(idxnan)
            continue
        end
        
        % demean & standardize
        tmpAlpha = (rawAlpha(j,~idxnan)-mean(rawAlpha(j,~idxnan)))./...
            std(rawAlpha(j,~idxnan));
        
        % return d/s alpha back to corresponding position in rawAlpha
        rawAlpha(j,~idxnan) = tmpAlpha;
    end
        
    % windsorization
    rawAlpha(rawAlpha > 3) = 3;
    rawAlpha(rawAlpha < -3) = -3;
end

% output
stdAlpha = rawAlpha;
end

%% ********************************************************************
function [trade,t0,lam,mu,booksize,tradesize] = myOptimizer(ymIdx,actStIdx,...
    sig,alphablend,volume,tcost,R,F)
% =====================================================================
% This function is used to perform optimization in Lecture 7
% Input:
%       sig - from myRisk()
%       alphablend - from myAlpha...()
%       Other parameter - either from or calculate from dataset PS2
% =====================================================================

booklimit = 100000000;  % 50*50 M book limit
tradelimit = 15000000;  % 15 M trade limit
[T,n] = size(tcost);
% get book size of each day in optimization that is used to select mu
booksize = size(T,3);   
% get trade size of each day in optimization that is used to select lamda
tradesize = size(T,3);

%setup options for the quadprog algorithm
options = optimset('Algorithm','interior-point-convex');
options = optimset(options,'Display','iter');

% give some value to lambda and mu
% lamda controls trade size, lager lamda smaller trade size
lam = 7;  
% mu controls book size, larger mu smaller book
mu = 0.004;       
rstar = 0.003*booklimit;  % 300,000 industry transaction limit per day for 50*50 M
fstar = 0.001*booklimit;  % 100,000 country transaction limit per day for 50*50 M
stTradeCap = 0.01*tradelimit;   % 150,000 sigle stock trade cap
avevolume = zeros(1,n);   % average daily volume
for i = 1:n
    idxnan = isnan(volume(:,i));
    avevolume(i) = mean(volume(~idxnan,i));
end

t0 = ymIdx(13); % first trading day
w = zeros(T,n); % initial position
trade = zeros(T,n);  % trading transaction

for i = t0:T        % daily time loop 
    disp(i)
    % at the beginning of the month, retrieve values for that month
    if ismember(i,ymIdx)
        % get month index
        [~,idxm] = ismember(i,ymIdx);
        % retrieve active stock index at the beginning of that month
        actIdx = actStIdx{idxm};
        if isempty(actIdx) 
            continue; 
        end
        % retrieve shrinkage estimator at the beginning of that month
        sigi = sig{idxm-12};
        % calculate beta from shrinkage estimator for that month
        % market index is the equal weighted composite of active stocks
        varm = sum(sum(sigi))./size(sigi,1)^2; % variance of market index
        covim = sum(sigi,2)./size(sigi,1);     % (i,m) covariance vector
        beta = covim./varm;      % beta vector
    end
    
    % retrieve alphablend for active stocks
    ai = alphablend(i,actIdx)';  
    % retrieve tcost for active stocks
    tcosti = tcost(i,actIdx)'; 
    
    % deactive stock if alpha or tcost is unavailable
    if any(isnan(ai)) || any(isnan(tcosti))
        idxnan = any([isnan(ai),isnan(tcosti)],2);
        actIdx(idxnan) = [];
        ai(idxnan) = [];
        tcosti(idxnan) = [];
        sigi = sigi(~idxnan,~idxnan);
        beta(idxnan) = [];
    end
        
    % build parameter matrix in daily 
    nact = size(actIdx,2);  % number of active stocks
    Ri = R(actIdx,:);       % industry dummy for active stocks
    Fi = F(actIdx,:);       % country dummy for active stocks
    wi = w(i-1,actIdx)';    % initial portfolio weight for active stocks
    
    % constraints for max trade size and position size (P13-14,Lecture 7)
    the = avevolume(actIdx)'.*0.02;  % trade size limit
    the(the > stTradeCap) = stTradeCap;
    pi = min(10.*the,0.025*booklimit/2);   % book size limit
    
    % matrix for objective formula (P26,Lecture 7)
    H = 2*mu.*[sigi,-sigi; -sigi,sigi];
    g = [2*mu.*sigi*wi-ai+lam.*tcosti; 
        -2*mu.*sigi*wi+ai+lam.*tcosti];
    
    % matrix for constrain (P27-29,Lecture 7)
    A = [Ri',-Ri'; -Ri',Ri'; Fi',-Fi'; -Fi',Fi']; 
    b = [rstar-Ri'*wi; rstar+Ri'*wi; fstar-Fi'*wi; fstar+Fi'*wi];
    C = [beta',-beta'];
    d = -beta'*wi;
    LB = zeros(2*nact,1);
    % relax trading constrain for the first day
    if i == t0
        UB = [max(0,pi-wi); max(0,pi+wi)];
    else
        UB = [max(0,min(the,pi-wi)); max(0,min(the,pi+wi))];
    end
    
    % run quadratic optimizer (P30,Lecture 7)
    [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    % if optimizer encounter error (-2), relax trade size constrain
    if exitflag == -2
        UB = [max(0,pi-wi); max(0,pi+wi)];
        [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    end
    
    % decompose u to y and z to get x and w (P22,Lecture 7)
    y = u(1:nact);  z = u(nact+1:end);  % long and short transaction
    x = wi+y-z;     % target weight
    % the target weight (x) will be the weight in the end of date i
    w(i,actIdx) = x;    
    trade(i,:) = w(i,:)-w(i-1,:);
    
    % book size and trade size of each day 
    booksize(i,1) = sum(abs(w(i,:)));   % total book size
    booksize(i,2) = sum(w(i,w(i,:) > 0));    % long book size
    booksize(i,3) = sum(w(i,w(i,:) < 0));    % short book size
    tradesize(i,1) = sum(abs(trade(i,:)));   % total trading size
    tradesize(i,2) = sum(trade(i,trade(i,:) > 0)); % long trading size
    tradesize(i,3) = sum(trade(i,trade(i,:) < 0)); % short trading size
end

end

function [back_weight,pnl,booksize,tradesize] = myBacktest(trade,ret,...
    tcost,myday,t0)
% =====================================================================
% This function is used to perform backtest by using trading transaction
% from optimization
% Input:
%       trade, t0 - from myOptimizer()
%       ret, tcost, myday - from given dataset in PS2 
% =====================================================================

[T,n] = size(trade);

% set NaN to 0 so that it will not affect calculation
trade(isnan(trade)) = 0;
ret(isnan(ret)) = 0;
tcost(isnan(tcost)) = 0;

back_weight = zeros(T,n);   % backtest weight of each day
pnld = zeros(T,1);  % PnL of each day
for i = t0:T
    back_weight(i,:) = back_weight(i-1,:).*(1+ret(i,:))+trade(i,:);
    pnld(i) = sum(back_weight(i-1,:).*ret(i,:))-sum(abs(trade(i,:)).*tcost(i,:));
end
pnl = cumsum(pnld); % cumulative PnL
booksize = sum(abs(back_weight),2); % book size in backtest of each day
tradesize = sum(abs(trade),2);      % trade size in backtest of each day

% build fints of cumulative PnL for plot
pnlts = fints(cellstr(myday(t0:end,:)),pnl(t0:end),'PnL','d');

% plot 
plot(pnlts)
set(gca,'FontSize',14)
title('Backtest from 02-Jan-1998 to 30-Dec-2002')
xlabel('Date')
ylabel('Cumulative PnL')
legend('hide')
end

function [sharpe,longest_dd,deepest_dd] = myStat(pnl,booksize,t0)
% =====================================================================
% This function is used to calculate statistics of backtest
% Input:
%       pnl, booksize - from myBacktest()
%       t0 - from myOptimizer()
% =====================================================================

T = size(pnl,1);
% daily return
ret = (pnl(t0+1:end)-pnl(t0:end-1))./booksize(t0:end-1);
retAnn = mean(ret)*252;      % annualized return
sigAnn = std(ret)*sqrt(252); % annualized volatility
sharpe = retAnn/sigAnn;     % annualized Sharpe ratio

wm = pnl(t0+1); % initialize water mark
ddt = 0;        % initialize drawdown time
dds = 0;        % initialize drawdown size
longest_dd = 0; % initialize longest drawdown time
deepest_dd = 0; % initialize longest drawdown size

for i = t0+2:T
    if pnl(i) > wm
        % increase water mark and set drawdown count to 0
        wm = pnl(i); ddt = 0; dds = 0;
    else
        % count drawdown
        ddt = ddt+1; dds = wm-pnl(i);
    end
    
    % compare with max drawdown time and size
    if ddt > longest_dd
        longest_dd = ddt;
    end
    if dds > deepest_dd
        deepest_dd = dds;
    end
end 
end
