function [ new_like ] = gmm_one_sample_posterior_like( gmm, x )
    [~, D] = size(x);
    
    % The MVNW part
    mu = (gmm.prior_scale * gmm.prior_mean + x) ./ (gmm.prior_scale + 1);
    dm = gmm.prior_mean - x;
    C  = gmm.prior_cov + (gmm.prior_scale / (gmm.prior_scale + 1))*(dm * dm');
    dof = gmm.prior_dof + 1;
    scale = gmm.prior_scale + 1;
    
    pred_dof = dof - D + 1;
    pred_cov = (scale + 1) / (scale * pred_dof) * C;
    
    mvtparams = make_mvt(mu, pred_cov, pred_dof);
    new_like  = fast_mvtpdf(mvtparams, x);

end

