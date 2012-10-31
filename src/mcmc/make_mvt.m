function [ mvtparams ] = make_mvt( mean_, S, dof )
%MAKE_MVT mean_ is N x D
    [N, D] = size(mean_);
    assert(N == 1);
    
    [mvtparams.istdev, t] = cholcov(inv(S));        
    assert(t == 0, 'S was not psd!');
    mvtparams.mean = mean_;
    mvtparams.dof = dof;

    mvtparams.pow = -(dof + D) / 2;
    
    % See the MATLAB documentation for Multivariate t Distribution
%     s = det(corr)^(-1/2);    
%     d = (dof * pi)^(-D/2);
%     g = gamma((dof + D) / 2) / gamma(dof / 2);

    log_s = log(det(mvtparams.istdev));
    assert(abs(log_s - log(sqrt(det(inv(S))))) < 1e-8);
    log_d = -D/2 * log(dof * pi);
    log_g = gammaln((dof + D) / 2) - gammaln(dof / 2);
    
%     mvtparams.Z = s * d * g;
    mvtparams.logZ = log_s + log_d + log_g;
%     assert(abs(mvtparams.Z - exp(mvtparams.logZ)) < 1e-8);

end

