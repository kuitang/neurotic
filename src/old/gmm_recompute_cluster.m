function [ gmm ] = gmm_recompute_cluster( gmm, X, k )
% Compute posterior predictive likelihood for idxs.
% Entire dataset is needed to compute parameters, but likelihoods are only
% computed for idxs.
% Uses gmm.K, gmm.s_z and gmm.n.
    assert(k > 0, 'k must identify a cluster!');
    if k == 1 % don't do anything for the background cluster
        return
    end
    [N, D] = size(X);
    kidxs  = gmm.s_z == k;
    Nk     = sum(kidxs);
    assert(Nk > 0, 'not computing an empty cluster!');    

    sample_mean = mean(X(kidxs,:), 1);
    sample_cov  = cov(X(kidxs,:), 1);
    sum_squares = Nk * sample_cov;
    
    gmm.scale(k)  = gmm.prior_scale + gmm.n(k);
    gmm.mean(k,:) = gmm.prior_scale .* gmm.prior_mean + Nk * sample_mean;
    gmm.mean(k,:) = gmm.mean(k,:) / gmm.scale(k);
    
    gmm.dof(k) = gmm.prior_dof + gmm.n(k);
    
    %% LAMBDA
    dm = gmm.prior_mean - sample_mean;
    km = gmm.prior_scale * gmm.n(k) / (gmm.prior_scale + gmm.n(k));
    gmm.cov(:,:,k) = gmm.prior_cov + sum_squares + km * (dm' * dm);
    [~, p] = cholcov(gmm.cov(:,:,k));
    assert(p == 0, 'gmm.cov(:,:,k) was not PD!');
    
    %% Posterior predictive
    % Murphy (228)
    gmm.pred_dof(k) = gmm.dof(k) - D + 1;
    gmm.pred_cov(:,:,k) = (gmm.scale(k) + 1) ./ (gmm.scale(k) .* gmm.pred_dof(k)) * gmm.cov(:,:,k);                

    gmm.pred_mvtparams{k} = make_mvt(gmm.mean(k,:), gmm.pred_cov(:,:,k), gmm.pred_dof(k));
    gmm.pred_x_like(:,k)  = fast_mvtpdf(gmm.pred_mvtparams{k}, X);
    assert(all(gmm.pred_x_like(:,k) > 0));
end

