function [ p ] = fast_mvtpdf( X, mvtparams )
    % We scale by both mean and covariance, which MATLAB
    % does not.    
    dX = bsxfun(@minus, X, mvtparams.mean);        
    Z = linsolve(mvtparams.cov_chol, dX', mvtparams.opts);
    
    log_p = mvtparams.logZ + mvtparams.pow * log((1 + sum(Z .* Z, 1) / mvtparams.dof))';
    
    p = exp(log_p);
end
