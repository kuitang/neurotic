function [ MK, SK ] = cond_moments( X, z )
% conditional MK(k,:) and SK(k,:,:) are conditional sample mean/cov
    [N, D] = size(X);
    K = max(z);
    
    MK = zeros(K, D);
    SK = zeros(D, D, K);
    for k = 1:K
        assert(sum(z==k) > 0);
        sX = X(z==k, :);
        MK(k,:)   = mean(sX, 1);
        SK(:,:,k) = cov(sX, 1); % cov(sX, 1) normalizes by N                
    end        

end

