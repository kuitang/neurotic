function [ h ] = sample_gmm_prior( gmm, K )


    X = [];    
    z = [];
    g = figure;
    for k = 1:K        
        C  = iwishrnd(gmm.prior_cov, gmm.prior_dof);
        mu = mvnrnd(gmm.prior_mean, C);
        myX = mvnrnd(mu, C, 100);
        subplot(1,K,k);
        hist(myX(:,3), 20);
        X = [ X ; myX ];        
        z = [ z ; k * ones(100, 1) ];
    end
    set(g, 'units', 'Normalized', 'position', [0 0.7 1 0.3]);
    h = figure;
    gscatter(X(:,1), X(:,2), z);
    set(h, 'units', 'pixels', 'position', [0 0 300 300]);
    xlim([0 200])
    ylim([0 200])

end

