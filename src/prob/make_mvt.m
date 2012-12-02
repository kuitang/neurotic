function [ mvtparams ] = make_mvt( mvtparams, mean_, cov_chol, dof )
%MAKE_MVT mean_ is N x D
    D = length(mean_);        
    assert(dof > 0);
        
    mvtparams.opts.UT = true;                
    mvtparams.cov_chol = cov_chol;    
    mvtparams.mean = mean_;
    mvtparams.dof = dof;
    mvtparams.pow = -(dof + D) / 2;
    
    % See the MATLAB documentation for Multivariate t Distribution
%     s = det(corr)^(-1/2);    
%     d = (dof * pi)^(-D/2);
%     g = gamma((dof + D) / 2) / gamma(dof / 2);
    
    % Cholesky decomposition fact:
    % log(det(cov)) == 2 * log(det(cov_chol)) == 2*sum*log(diag(cov_chol))
    log_s = -sum(log(diag(cov_chol)));
    
%    assert(abs(log_s - log(sqrt(det(inv(S))))) < 1e-8);
    log_d = -D/2 * log(dof * pi);
    log_g = gammaln((dof + D) / 2) - gammaln(dof / 2);
    
%     mvtparams.Z = s * d * g;
    mvtparams.logZ = log_s + log_d + log_g;
%     assert(abs(mvtparams.Z - exp(mvtparams.logZ)) < 1e-8);

end

