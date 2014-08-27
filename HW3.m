clear;
load('C:\Users\maitreyi\Documents\MATLAB\database_v04.mat')


% %3a
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
F=zeros(n,length(uni_coun));
for i=1:n
    for j=1:length(uni_coun)
        F(i,j)=strcmp(uni_coun{j},countries{i});
    end
end


T=1504;

w=zeros(21,1);
for i=1:21
    w(22-i)=1/11-1/231*(i-1);
end
rt=zeros(T,n);
rt(2:T,:)=tri(2:T,:)./tri(1:T-1,:)-1;
alpharev_=zeros(T,n);
alpharev=zeros(T,n);

% for i=1:T
%     for j=1:n
%         if isnan(rt(i,j))
%             rt(i,j)=0;
%         end
%     end
% end
% 
% for i=22:T
%     alpharev_(i,:)=-w'*rt(i-20:i,:);
% end
% 
RRRR=R*inv(R'*R)*R';
% for i=22:T
%     alpharev(i,:)=alpharev_(i,:)*(eye(n)-RRRR);
% end
% 
% %%3b
% alpharec=zeros(T,n);
% w=zeros(45,1);
% for i=1:45
%     w(46-i)=1/23-1/1035*(i-1);
% end
% for i=46:T
%     alpharec(i,:)=-w'*rec(i-44:i,:);
% end
% 
% %%3c
% alphaval=mtbv;
% %3d
% alphamom=zeros(T,n);
% alphamom(232:T,:)=tri(232:T,:)./tri(1:T-231,:)-1;
% 
% ma=mean(alpharev,2);
% mb=mean(alpharec,2);
% mc=mean(alphaval,2);
% md=mean(alphamom,2);
% stda=std(alpharev,0,2);
% stdb=std(alpharec,0,2);
% stdc=std(alphaval,0,2);
% stdd=std(alphamom,0,2);
% 
% for i=1:n
%     alpharev(:,i)=(alpharev(:,i)-ma)./stda;
%     alpharec(:,i)=(alpharec(:,i)-mb)./stdb;
%     alphaval(:,i)=(alphaval(:,i)-mc)./stdc;
%     alphamom(:,i)=(alphaval(:,i)-md)./stdd;
% end
% 
% alpharev=max(alpharev,-3);
% alpharev=min(alpharev,3);
% alpharec=max(alpharec,-3);
% alpharec=min(alpharec,3);
% alphaval=max(alphaval,-3);
% alphaval=min(alphaval,3);
% alphamom=max(alphamom,-3);
% alphamom=min(alphamom,3);
% 
% alphablend=.5*alpharev+.25*alpharec+.15*alphaval+.10*alphamom;
% me=mean(alphablend,2);
% stde=std(alphablend,0,2);
% for i=1:n
%     alphablend(:,i)=(alphablend(:,1)-me)./stde;
% end
% alphablend=max(alphablend,-3);
% alphablend=min(alphablend,3);
% 
% %
% %4
% load shrink.mat
% myday=cellstr(myday);
% monthposition=zeros(72,1);
% monthposition(1)=1;
% pos=2;
% prem=1;
% for i=1:T
%     if month(myday(i))~=prem
%         monthposition(pos)=i;
%         pos=pos+1;
%     end
%     prem=month(myday(i));
% end
% 
% for i=1:60
%     S=cov(rt(monthposition(i):monthposition(i+12),:));
%     sigma_=mean(var(rt(monthposition(i):monthposition(i+12),:)));
%     SIGMA{i}=(1-shrink(i))*sigma_*eye(n)+shrink(i)*S;
% end
% 
% 
% for i=1:60
%     Beta{i}=zeros(n,1);
%     mktrtn=mean(rt(monthposition(i):monthposition(i+12),:),2);
%     for j=1:n
    %     Beta{i}(j)=regress(rt(monthposition(i):monthposition(i+12),j),mktrtn);
%     end
% end
% 
lambda=2;
mu=0.002;
% pi
% theta
options=optimset('Algorithm','interior-point-convex');
options=optimset(options,'Display','off');
r=300000;
f=100000;
% A=[R' -R'; -R' R'; F' -F'; -F' F'];
% LB=zeros(2*n,1);
% w=zeros(n,1);
% x=zeros(n,1);
% C=[-Beta{1} Beta];
% d=-Beta{1}'*w;
% H=2*mu*[SIGMA{1}, SIGMA{1}; SIGMA{1}, SIGMA{1}];
% g=[2*mu*SIGMA{1}*w-alphablend(1,:)'+lambda*tcost(1,:)'; ...
%     -2*mu*SIGMA{1}*w+alphablend(1,:)'+lambda*tcost(1,:)'];
% b=[r-R'*w; r+R'*w; f-F'*w; f+F'*w];


%         act_tri = tri(max(1,t-253):t-1,:);
%         act_tri_nan = sum(isnan(act_tri));
%         act_tri_nan = find(act_tri_nan~=0);
%         if(~isempty(act_tri_nan))
%             active_list(1,act_tri_nan) = 0;
%         end
%         %
%         act_volume = volume(max(1,t-253):t-1,:);
%         act_volume_nan = sum(isnan(act_volume));
%         act_volume_nan = find(act_volume_nan~=0);
%         if(~isempty(act_volume_nan))
%             active_list(1,act_volume_nan) = 0;
%         end
%         %
%         act_tcost = tcost(max(1,t-253):t-1,:);
%         act_tcost_nan = sum(isnan(act_tcost));
%         act_tcost_nan = find(act_tcost_nan~=0);
%         if(~isempty(act_tcost_nan))
%             active_list(1,act_tcost_nan) = 0;
%         end
%         %
%         act_ret = tri(max(1,t-254)+1:t-1,:)./tri(max(1,t-254):t-2,:) - 1; 
%         act_ret_nan = sum(isnan(act_ret));
%         act_ret_nan = find(act_ret_nan~=0);
%         if(~isempty(act_ret_nan))
%             active_list(1,act_ret_nan) = 0;
%         end
%         %




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
trade_days = int16(size(price,1));
stk_num = int16(size(allstocks,2));
book_size = zeros(trade_days,1);
trade_size = zeros(trade_days,1);
back_weight = zeros(trade_days,stk_num);
pnl = zeros(trade_days,stk_num);
trade_time = int16(datevec(myday));
trade_time(:,4:6) = [];
start_trade_idx = int16(find(trade_time(:,1)==1998 & trade_time(:,2)==1 & trade_time(:,3)==2));
cur_month = int16(trade_time(start_trade_idx,2));
last_month = int16(trade_time(start_trade_idx,2));
active_list = isactivenow(start_trade_idx,:);
t = start_trade_idx;
act_tri = tri(max(1,t-253):t-1,:);
act_tri = act_tri(:,active_list);
act_industry = R(active_list',:);
act_rec = rec(max(1,t-253):t-1,:);
act_rec = act_rec(:,active_list);
act_ret = tri(max(1,t-254)+1:t-1,:)./tri(max(1,t-254):t-2,:) - 1;
act_tcost = tcost(max(1,t-253):t-1,:);
act_tcost = act_tcost(:,active_list);
act_ret = act_ret(:,active_list);
act_R = R(active_list,:);
act_RRRR = RRRR(active_list,active_list);
act_F = F(active_list,:);
act_stk_num = sum(active_list);
act_intpos = back_weight(t-1,:);
act_intpos = act_intpos(:,active_list);
act_volume = volume(max(1,t-253):t-1,:);
act_volume = act_volume(:,active_list);
shrink = zeros(60,1);
alpharev = zeros(trade_days,stk_num);
alpharec = zeros(trade_days,stk_num);
alphaval = zeros(trade_days,stk_num);
alphamom = zeros(trade_days,stk_num);
alphablend = zeros(trade_days,stk_num);
t0 = start_trade_idx;
trade = zeros(trade_days,stk_num);

% shrink_idx = 1;
% shrink_cov = SIGMA(shrink_idx);
% loop through every day
for t = start_trade_idx:trade_days
    % extract not NaN and active stocks
    cur_month = trade_time(t,2);
    if cur_month~=last_month
        last_month = cur_month;
%         shrink_idx = shrink_idx + 1;
%         shrink_cov = SIGMA(shrink_idx);
        %
        active_list = isactivenow(t,:);
    end   
    
    act_tri = tri(max(1,t-253):t-1,:);
    act_tri = act_tri(:,active_list);
    act_tri(isnan(act_tri)) = 0;
    act_volume = volume(max(1,t-253):t-1,:);
    act_volume = act_volume(:,active_list);
    act_volume(isnan(act_volume)) = 0;
    act_tcost = 1/10000 + 0.5*tcost(max(1,t-253):t-1,:);
    act_tcost = act_tcost(:,active_list);
    act_tcost(isnan(act_tcost)) = 0;
    act_ret = tri(max(1,t-254)+1:t-1,:)./tri(max(1,t-254):t-2,:) - 1;        
    act_ret = act_ret(:,active_list);
    act_ret(isnan(act_ret)) = 0;
    %
    act_rec = rec(max(1,t-253):t-1,:);
    act_rec = act_rec(:,active_list);        
    act_R = R(active_list,:);
    act_RRRR = RRRR(active_list,active_list);
    act_F = F(active_list,:);
    act_stk_num = sum(active_list);
    act_intpos = back_weight(t-1,:);
    act_intpos = act_intpos(:,active_list);
    act_industry = R(active_list',:);
    % calculate shrinkage estimator
    act_cov = shrinkage_estimator(act_ret);
%     act_cov = shrink_cov(active_list',:);
%     act_cov = shrink_cov(:,active_list);
    % calculate alpha
    w=zeros(21,1);
    for i=1:21
        w(22-i)=1/11-1/231*(i-1);
    end
    alpha_a = (-w'*act_ret(end-20:end,:))*(eye(act_stk_num)-act_RRRR);
    w=zeros(45,1);
    for i=1:45
        w(46-i)=1/23-1/1035*(i-1);
    end
    alpha_b = -w'*rec(t-44:t,:);
    alpha_b = alpha_b(:,active_list);
    alpha_c = mtbv(t, active_list);
    alpha_d = act_tri(1,:)./act_tri(end,:) - 1;
    % demean, standard and windsorization
    std_alpha_a = standardize(alpha_a,ones(size(alpha_a)),2);
    std_alpha_b = standardize(alpha_b,ones(size(alpha_b)),2);
    std_alpha_c = standardize(alpha_c,ones(size(alpha_c)),2);
    std_alpha_d = standardize(alpha_d,ones(size(alpha_d)),2);
    std_alpha = 0.5*std_alpha_a + 0.25*std_alpha_b + 0.15*std_alpha_c + 0.1*std_alpha_d;
    std_alpha = standardize(std_alpha,ones(size(std_alpha)),2);
    alpharev(t,active_list) = std_alpha_a;
    alpharec(t,active_list) = std_alpha_b;
    alphaval(t,active_list) = std_alpha_c;
    alphamom(t,active_list) = std_alpha_d;
    alphablend(t,active_list) = std_alpha;
    % calculate weights by running optimization
    b = [r-act_R'*act_intpos'; r+act_R'*act_intpos'; f-act_F'*act_intpos'; f+act_F'*act_intpos'];
    H = 2*mu*[act_cov,-act_cov;-act_cov,act_cov];
    g = [2*mu*act_cov*act_intpos'-std_alpha'+lambda*act_tcost(end,:)';
        -2*mu*act_cov*act_intpos'+std_alpha'+lambda*act_tcost(end,:)'];
    act_mkt_ret = mean(act_ret,2);
    beta = zeros(act_stk_num,1);
    for i = 1:act_stk_num
        tmp = regress(act_ret(:,i),[ones(size(act_mkt_ret)),act_mkt_ret]);
        beta(i,1) = tmp(2);        
    end
    C = [beta' -beta'];
    d = -beta'*act_intpos';
    A = [act_R',-act_R';-act_R',act_R';act_F',-act_F';-act_F',act_F'];
    LB = zeros(2*act_stk_num,1);
    % WTF
    act_volume(isnan(act_volume)) = 0;
    theta = min(0.01*act_volume(end,:),150000*ones(1,act_stk_num));
%     theta = min(0.01*volume(t-1,idx)', 150000*ones(n_active,1));
    pie = min(10*theta,1250000);
    gamma = max(act_intpos-theta,-pie);
    sigmma = min(act_intpos+theta,pie);
    UB = [max(0,min(theta,pie-act_intpos))';
        max(0,min(theta,pie+act_intpos))'];
    [u,fval,exitflag,output] = quadprog(H,g,A,b,C,d,LB,UB,[],options);
    u = transpose(u);
    if(isempty(u))
        x = zeros(1,act_stk_num);
    else
        x = - u(1,1:act_stk_num) + u(1,act_stk_num+1:end);    
    end
    % calculate profits    
    total_ret = tri(t,:)./tri(t-1,:)-1;
    total_ret(1,isnan(total_ret)) = 0;
    back_weight(t,:) = back_weight(t-1,:).*(1+total_ret);
    back_weight(t,active_list) = back_weight(t,active_list) + x;
    pnl(t,:) = back_weight(t-1,:).*total_ret;
    pnl(t,active_list) = pnl(t,active_list) - abs(x.*act_tcost(end,:));
    fprintf('%s %d %d %d\n','Hello World!',size(u,1),size(u,2),sum(active_list));
end