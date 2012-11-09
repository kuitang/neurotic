function [ x_like ] = gmm_pred_x_like( gmm, x )
    x_like = zeros(gmm.K, 1);
    x_like(1) = 1/gmm.N * gmm.background_pdf(x(3));
    
    for kk = 2:gmm.K
        k = gmm.k_idx(kk);
        x_like(k) = fast_mvtpdf(gmm.pred_mvtparams{k}, x);
    end

end

