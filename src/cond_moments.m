function [ M, S, NM, NS ] = cond_moments( X, K, z )
% conditional MK(k,:) and SK(k,:,:) are conditional sample mean/cov
    [Nk, D] = size(X);    
    
    M  = zeros(K, D);
    NM = zeros(K, D);
    S = zeros(D, D, K);
    NS = zeros(D, D, K);
    for k = 1:K
        Nk = sum(z==k);        
        if Nk > 0
            sX = X(z==k, :);
            M(k,:)    = mean(sX, 1);
            NM(k,:)   = Nk*M(k,:);
            S(:,:,k)  = cov(sX, 1); % cov(sX, 1) normalizes by N
            NS(:,:,k) = Nk*S(:,:,k);                    
        end
    end        

end

