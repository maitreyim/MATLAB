% optimization and backtest

clear all; clc;

load database_v04
load alpha


T = size(rec,1);
n = size(rec,2);
ret = zeros(T,n);
ret(2:end,:) = tri(2:end,:)./tri(1:end-1,:)-1;
t0 = 246;
% t0 = month_start(13);  % first day of trading: Jan 1, 1998
shrink = zeros(60,1);

% clear dataset
for t = t0-1:T
    nodata = isnan(tcost(t,:));
    tcost(t,nodata) = nanmean(tcost(t-23:t-1,nodata));
    nodata = isnan(alphablend(t,:));
    alphablend(t,nodata) = 0;
end


trade = zeros(T,n);
back_weight = zeros(T,n);
pnl = zeros(T,1);
booksize = zeros(T,1);
tradesize = zeros(T,1);


lambda = 4;
mu = 0.003;

trade_time = datevec(myday);
trade_time(:,4:6) = [];
start_trade_idx = find(trade_time(:,1)==1998 & trade_time(:,2)==1 & trade_time(:,3)==2);
cur_month = trade_time(start_trade_idx,2);
last_month = trade_time(start_trade_idx,2);

n=566;
for i=1:n
    industries{i}=allstocks(1,i).industrylist.industry;
end
industry=unique(industries);
rho=length(industry);
R=zeros(n,rho);
for i=1:n
    for j=1:rho
        R(i,j)=strcmp(industry{j},industries{i});
    end
end
industry = R;

countries = cell(n,1);
for i=1:n
    if(size(allstocks(1,i).indexlist,1)~=1)
        for j=1:size(allstocks(1,i).indexlist,1)
            if(~isempty(allstocks(1,i).indexlist(j).index))
                allstocks(1,i).indexlist = struct('date','01-Jan-1998','index',allstocks(1,i).indexlist(j).index);
                break;
            end

        end
    end
end
for i=1:n
    countries(i,1)=cellstr(allstocks(1,i).indexlist.index);
end
uni_coun = unique(countries);
country=zeros(n,length(uni_coun));
for i=1:n
    for j=1:length(uni_coun)
        country(i,j)=strcmp(uni_coun{j},countries{i});
    end
end
shrink_idx = 1;
idx = find(isactivenow(t0,:));  

for t = t0:T
    t;
    cur_month = trade_time(t,2);
    if cur_month~=last_month % if t is the first day of a month, reset the universe   
        idx = find(isactivenow(t,:));  
        last_month = cur_month;
        shrink(shrink_idx,1) = betaa;
        shrink_idx = shrink_idx + 1;
    end
    n_active = length(idx);

    % position as of close on day t-1
    w = back_weight(t-1,idx)';  % what if stocks are active in previous month but inactive this month?
    
    % alpha up to day t-1
    alpha = alphablend(t-1,idx)';
    
    % transaction cost model
    % 1bp+0.5*bid-ask spread
    tau = 1/10000+0.5*tcost(t-1,idx)';
    
    % risk model
    % shrinkage estimator of covariance matrix of stock returns
    offset = (year(myday(t,:))-1997)*12+month(myday(t,:))-12;
%     sigma = cov{offset};    
    act_ret = tri(max(1,t-254)+1:t-1,:)./tri(max(1,t-254):t-2,:) - 1;        
    act_ret = act_ret(:,idx);
    act_ret(isnan(act_ret)) = 0;
    [sigma,betaa] = shrinkage_estimator(act_ret);
    fprintf('%s\n','Hi');
    
    % constraints
    theta = min(0.01*volume(t-1,idx)', 150000*ones(n_active,1));   % maximum trade size
    pi = min(10*theta, 1250000*ones(n_active,1));                  % maximum position size
%     beta = betalist(offset,idx)';                                  % market exposure
    beta = random('normal',1,0.5,[n_active,1]);

%     beta = zeros(n_active,1);
%     for i = 1:act_stk_num
%         tmp = regress(act_ret(:,i),[ones(size(act_mkt_ret)),act_mkt_ret]);
%         beta(i,1) = tmp(2);        
%     end
    
    R = industry(idx',:);                                          % industry exposure
    F = country(idx',:);                                           % country exposure
    
    rstar = 300000;                                                % limit of industry exposure
    fstar = 100000;                                                % limit of country exposure
    
    % quadratic programming
    H  = 2*mu*[sigma -sigma; -sigma sigma];
    g  = [2*mu*sigma*w-alpha+lambda*tau; -2*mu*sigma*w+alpha+lambda*tau];
    A  = [R' -R'; -R' R'; F' -F'; -F' F'];
    b  = [rstar-R'*w; rstar+R'*w; fstar-F'*w; fstar+F'*w];
    C  = [beta' -beta'];
    d  = -beta'*w;
    LB = zeros(2*n_active,1);
    UB = [max(zeros(n_active,1),min(theta,pi-w)); max(zeros(n_active,1),min(theta,pi+w))];

    options = optimset('Algorithm','interior-point-convex');
    [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    while exitflag == -2
        theta = theta*2;
        pi = min(10*theta, 1250000*ones(n_active,1));
        UB = [max(zeros(n_active,1),min(theta,pi-w)); max(zeros(n_active,1),min(theta,pi+w))];
        [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    end
    
    fprintf('%d %d\n',size(u,1),size(u,2));
    
    y = u(1:n_active,1);
    z = u(n_active+1:end,1);
    trade(t,idx) = y-z;
    
    
    % backtest
    kickout = isactivenow(t-1,:) & (~isactivenow(t,:));                  % active stocks becoming inactive
    trade(t,:) = trade(t,:)-back_weight(t-1,:).*(1+ret(t,:)).*kickout;   % clear all positions in stocks that is kicked out
    
    back_weight(t,:) = back_weight(t-1,:).*(1+ret(t,:)) + trade(t,:);    % weight today after trade
    back_weight(t,isnan(back_weight(t,:))) = 0;                          % clear positions that is NaN due to w
    
    pnl(t) = pnl(t-1) + nansum(back_weight(t-1,:).*ret(t,:))-nansum(abs(trade(t,:)).*(1/10000+0.5*tcost(t,:)));   % cumulative p/l

    booksize(t) = nansum(abs(back_weight(t,:)));
    tradesize(t) = nansum(abs(trade(t,:)));

end

shrink(shrink_idx) = beta;

% Sharpe ratio
dpnl = diff(pnl);       % daily P/L
sharpe = mean(dpnl(t0:T-1))/std(dpnl(t0:T-1))*sqrt(252);

% longest drawdown
dd = zeros(T,1);
for i = t0:T
    for j = i:T
        if pnl(i) < pnl(j)
            dd(i) = j-i;
            break
        end
    end
end
longest_dd = max(dd);

% deepest drawdown
dd = zeros(T,1);
for i = t0:T
    dd(i) = max(pnl(t0:i))-pnl(i);
end
deepest_dd = max(dd);



save('backtest.mat','shrink','alpharev','alpharec','alphaval','alphamom',...
    'alphablend','lambda','mu','t0','trade','back_weight','pnl','booksize',...
    'tradesize','sharpe','longest_dd','deepest_dd')
