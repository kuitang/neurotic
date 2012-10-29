function [ MX ] = cond_sum( X, z )
%cond_sum conditional sum by responsibility z: SX(k) = sum(X(z==k,:))
    [N, D] = size(X);
    K = max(z);
    
    MX = zeros(K, D);
    for k = 1:K
        idxs = z==k;
        if sum(idxs) == 0
          MX(k,:) = zeros(1, D);
        else
          MX(k,:) = mean(X(idxs, :), 1);
        end
    end        

end

