% =====================================================================
% This code is created by Group 5: Jie (Edward) Sheng, Anshul Maheshwari
% Chao Liang, Li Zhang for MGMT237M1 Problem Set 2
% =====================================================================

% initialize working environment
clear; close all; clc;
% add new working path if necessary
path(path, '')

% Step 1. import allstocks provided by professor
load('allstocks.mat');  % it contains allstocks structure from PS1
n = size(allstocks,2);  % n is number of unique dscode
% asdscode is the list of dscode we are looking for. It will be used
% to extract infomation from TickSummary.mat
asdscode = cell(size(allstocks));
for i = 1:size(asdscode,2)
    asdscode{i} = allstocks(i).dscode;
end
xlswrite('MGMT237M1_HW2.xlsm',asdscode,'asdscode')

% asibes is the list of IBES Ticker in allstocks that will be used to 
% extract recommendation from WRDS
asibes = cell(size(allstocks));
for i = 1:size(asibes,2)
    asibes{i} = allstocks(i).ibeslist(1).ibes;
end
xlswrite('MGMT237M1_HW2.xlsm',asibes,'asibes')

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 2. import data from DataStream
% Reformatted Data are store in MGMT237M1_HW2.xlsm. Date series are 
% reformatted by macros and are ready to be input into MATLAB.  
% Detail description can be find in Note worksheet in MGMT237M1_HW2.xlsm or
% in written document. 
% Data are stored in the same order of asdscode as column and myday as row
% unadjusted price in local currency
% create myday (T×1) variable
[~,myday] = xlsread('MGMT237M1_HW2.xlsm','myday');
T = size(myday,1);

price = xlsread('MGMT237M1_HW2.xlsm','UP'); 
% total return index in Euros
tri = xlsread('MGMT237M1_HW2.xlsm','XRIE');
% unadjusted volume in 1,000 shares
UVO = xlsread('MGMT237M1_HW2.xlsm','UVO');
% unadjusted price in Euros
XUPE = xlsread('MGMT237M1_HW2.xlsm','XUPE');
% market value to book value
mtbv = xlsread('MGMT237M1_HW2.xlsm','MTBV');
% market capitalization in million Euros
cap = xlsread('MGMT237M1_HW2.xlsm','XMVE');

% adjust data to meet the requirement
% daily volume in Euros
volume = UVO.*XUPE.*1000;
% market capitalization in Euros
cap = cap.*1000000;

% create asmon which are the number list of the year and month of myday
% that will be used to compare with monthlist in TickSummary
asmon = xlsread('MGMT237M1_HW2.xlsm','asmon');

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 3. import TickSummary provided by professor
load('TickSummary.mat'); 
% this file contains three variables: dscode, monthlist, and spread
% index dscode from allstocks to dscode from TickSummary to extract data
% from TickSummary
[~,idxds] = ismember(asdscode,dscode);
% index date
monthlist(:,3) = []; % only use year and month column of monthlist
[~,idxday] = ismember(asmon,monthlist,'rows');
% extract bid-ask spread from TickSummary to build bidask matrix
tcost = NaN(T,n);
for i = 1:T     % i controls time loop
    for j = 1:n     % j controls stock loop
        if all([idxday(i),idxds(j)])    % make sure index is valid
            % take average of all 7 b-a spreads for each stock on each day
            tcost(i,j) = mean(spread(idxday(i),idxds(j),1:end));
        end
    end
end

% backfill missing data with the first available data
for j = 1:n     % j controls stock loop
    idx = find(~isnan(tcost(:,j)),1,'first');
    tcost(1:idx-1,j) = tcost(idx,j);
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 4. build recomendation matrix
% read recommendation raw data, raw data has AggUpgrade which is total 
% upgrade - total downgrade.  It has duplicate date that need to be deleted
% and 0 sum that need to be deleted
recraw = dataset('XLSFile','MGMT237M1_HW2.xlsm','Sheet','rectext');
% delete duplicate date that has multiple signal at that day, the last day
% of that duplicate date has AggUpgrade that is the total upgrade - total
% downgrade
recraw = unique(recraw,{'IBES','Datetext'},'last');
% delete 0 sum
recraw(recraw.AggUpgrade == 0,:) = [];

% noninal IBES for comparison
recraw.IBES = nominal(recraw.IBES);

% build rec matrix
rec = zeros(T,n);
for j = 1:n     % j controls stock loop
    idxrec = recraw(recraw.IBES == asibes{j},:);    % screen rec for this stock
    % get index of rec in daily time series
    [~,idxrec.idx] = ismember(idxrec.Datetext,myday);
    % remove invalid idx
    idxrec(idxrec.idx == 0,:) = [];
    % extract AggUpgrade to corresponding place in rec matrix
    rec(idxrec.idx,j) = idxrec.AggUpgrade;
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 5. construct isactivenow matrix
isactivenow = zeros(T,n);

% get previous month that can be used for index previous month for criteria e)
[~,asmontext] = xlsread('MGMT237M1_HW2.xlsm','asmontext');
currmon = asmontext(:,1);
premon = asmontext(:,2);
[~,idxpremon] = ismember(premon,currmon);

for j = 1:n     % j controls stock loop
    for i = 253:T   % i controls time loop
        criteria = zeros(1,6);  % criteria matrix
        
        % criteria a)
        if sum(isnan(price(i-252:i-1,j)))/252 < 0.1 && ...
                sum(isnan(tri(i-252:i-1,j)))/252 < 0.1
            criteria(1) = 1;
        end
        
        % criteria b)
        if ~any(isnan(price(i-10:i-1,j)))
            if range(price(i-10:i-1,j)) 
                criteria(2) = 1;
            end
        end
        
        % criteria c)
        if sum(isnan(volume(i-21:i-1,j)))/21 < 0.1 && ...
                sum(isnan(tri(i-21:i-1,j)))/21 < 0.1
            criteria(3) = 1;
        end
        
        % criteria d)
        if ~isnan(mtbv(i,j)) && ~isnan(cap(i,j))
            criteria(4) = 1;
        end
        
        % criteria e)
        if idxpremon(i)     % validate index
            if ~isnan(tcost(idxpremon(i),j))
                criteria(5) = 1;
            end
        end
        
        % criteria f)
        if range(rec(1:i-1,j))
            criteria(6) = 1;
        end
        
        % result
        if all(criteria)
            isactivenow(i,j) = 1;
        end
    end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Step 6. refine isactivenow
% find the index of the beginning of the month
[~,idxmonstart] = unique(currmon,'first');
idxmonstart = sort(idxmonstart);
% find the index of the end of the month
[~,idxmonend] = unique(currmon);
idxmonend = sort(idxmonend);

for j = 1:n     % j controls stock loop
    for i = 2:size(idxmonstart,1)   % i controls time loop
        % index start/end of this month and last month
        monstart = idxmonstart(i); premonstart = idxmonstart(i-1);
        monend = idxmonend(i); premonend = idxmonend(i-1);
    
        filter = zeros(1,2);
        
        % first filter
        if (isactivenow(monstart-1,j) && tcost(monstart,j) < 0.01) || ...
                (~isactivenow(monstart-1,j) && tcost(monstart,j) < 0.008)
            filter(1) = 1;
        end

        % second filter
        if (isactivenow(monstart-1,j) && ...
                mean(volume(premonstart:premonend,j)) > 1000000) || ...
                (~isactivenow(monstart-1,j) && ...
                mean(volume(premonstart:premonend,j)) > 1200000)
            filter(2) = 1;
        end
        
        % if do not pass both filter, inactive stock for the whole month
        if ~all(filter)
            isactivenow(monstart:monend,j) = 0;
        end
    end
end

whos allstocks myday price tri volume mtbv cap tcost rec isactivenow
save('MGMT237M1_HW2_Group5','allstocks','myday','price','tri','volume',...
    'mtbv','cap','tcost','rec','isactivenow')
