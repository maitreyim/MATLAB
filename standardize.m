% std_alpha = standardize(raw_alpha,isactivenow,iter_num)
% input: 
%   raw_alpha, [T,n] matrix
%   isactivenow, [T,n] matrix; iter_num, scalar
% output: 
%   std_alpha, [T,n] matrix

function std_alpha = standardize(raw_alpha,isactivenow,iter_num)
    % initialization
    if size(raw_alpha)~=size(isactivenow)
        fprintf('%s\t','The input dimension do not match');
        exit;
    end
    std_alpha = raw_alpha;
    [~,col_num] = size(isactivenow);
    for i = 1:iter_num
        % demean
        average = repmat(sum(std_alpha,2)./sum(isactivenow,2),1,col_num);
        std_alpha = std_alpha - average.*isactivenow;
        % standardize variance
        variance = repmat(sqrt(sum(std_alpha.^2,2)./(sum(isactivenow,2)-1)),1,col_num);
        std_alpha = std_alpha./variance;
        % windsorization
        std_alpha(std_alpha>3) = 3;
        std_alpha(std_alpha<-3) = -3;
    end
end