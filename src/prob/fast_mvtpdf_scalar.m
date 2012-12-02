function [ p ] = fast_mvtpdf_scalar( x, mvtparams )
    dX = x - mvtparams.mean;
    Z = linsolve(mvtparams.cov_chol, dX', mvtparams.opts);
    
    log_p = mvtparams.logZ + mvtparams.pow * log((1 + sum(Z .* Z, 1) / mvtparams.dof))';
    
    p = exp(log_p);
end
