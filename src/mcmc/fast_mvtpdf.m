function [ p ] = fast_mvtpdf( mvtparams, X )
    % We scale by both mean and covariance, which MATLAB
    % does not.    
    dX = bsxfun(@minus, X, mvtparams.mean);
    dX = dX * mvtparams.istdev'; % MATLAB's chol is L'L!
    
    %p = mvtparams.Z * (1 + sum(dx.^2) / mvtparams.dof) ^ mvtparams.pow;
    log_p = mvtparams.logZ + mvtparams.pow * log((1 + sum(dX .* dX, 2) / mvtparams.dof));
    %assert(abs(p - exp(log_p)) < 10e-8);
    %assert(p > 0);
    p = exp(log_p);
end
