% shrinkage covariance matrix estimator
% input:
%   data, [n,p] matrix, p random variables and n observations
% output:
%   shrinkage_estimator, [p,p] matrix
%   shrinkage_intensity, (1-beta)*m_n

function [shrinkage_estimator, shrinkage_intensity] = shrinkage_estimator(data)
    % initializtion
    [row_num, col_num] = size(data);
    % naive estimator
    naive_estimator = cov(data);
    identity = eye(col_num);
    % estimate \mu_n by m_n
    m_n = trace(naive_estimator*identity)/col_num;
    % estimte \beta^2
    beta = 0;
    for i = 1:row_num
        beta = beta + sum(sum((data(i,:)'*data(i,:)-naive_estimator).^2));
    end
    beta = beta/sum(sum((naive_estimator-m_n*identity).^2));
    beta = 1 - beta/((row_num)*(row_num-1));
    % output    
    shrinkage_intensity = beta;
    shrinkage_estimator = shrinkage_intensity*identity + beta*naive_estimator;
end