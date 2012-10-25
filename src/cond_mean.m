function [ MX ] = cond_sum( X, z )
%cond_sum conditional sum by responsibility z: SX(k) = sum(X(z==k,:))
    [N, D] = size(X);
    K = max(z);
    
    MX = zeros(K, D);
    for k = 1:K
        assert(sum(z==k) > 0);
        MX(k,:) = mean(X(z==k, :), 1);
    end        

end

