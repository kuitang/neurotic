function [ mvtparams ] = make_mvt( mean_, corr, dof )
%MAKE_MVT Same as MATLAB params, with scaled mean and corr.
    [N, D] = size(mean_);
    assert(N == 1);
    
    [stdev, neg_eigs] = cholcov(corr);
    assert(neg_eigs == 0, 'corr was not positive definite!');
    
    mvtparams.mean = mean_;
    mvtparams.dof = dof;
    mvtparams.stdev = chol(inv(corr));
    mvtparams.pow = -(dof + D) / 2;
    
    % See the MATLAB documentation for Multivariate t Distribution
    s = det(corr)^(-1/2);    
    d = (dof * pi)^(-D/2);
    g = gamma((dof + D) / 2) / gamma(dof / 2);
    
    mvtparams.Z = s * d * g;
    

end

