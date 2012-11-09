function [ h ] = sample_gmm_prior( gmm, K )


    X = [];    
    z = [];
    for k = 1:K
        C  = iwishrnd(gmm.prior_cov, gmm.prior_dof);
        mu = mvnrnd(gmm.prior_mean, C);
        X = [ X ; mvnrnd(mu, C, 100) ];
        z = [ z ; k * ones(100, 1) ];
    end
    h = figure;
    gscatter(X(:,1), X(:,2), z);
    set(h, 'position', [0 0 300 300], 'units', 'pixels');
    xlim([0 200])
    ylim([0 200])

end

