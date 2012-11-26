function [ samples ] = gmm_gibbs_diag(X, img, gmm, niters, burn_in, ...
                                      sample_freq, recovery_file )
% gmm_gibbs_diag(X, img, gmm, niters, burn_in, sample_freq, recovery_file )
    assert(nargin >= 5, 'not enough arguments!');
    
    [N, D] = size(X);
    
    samples = [];

    if nargin >=7 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end        
    
    figure
    
    tic
    for i = first:niters
        tic
        z_old = gmm.z;
        gmm = gmm_gibbs_diag_iter(gmm, X);
        change = mean(abs(z_old ~= gmm.z));
        disp(['Iteration ' num2str(i) ' took ' num2str(toc) ... 
              ' seconds; ' num2str(change) ' changed classes']);        

        if i > burn_in && mod(i, sample_freq) == 0
            samples = [samples gmm];
            iter = i;
            save(['iter' num2str(i) '.mat']);
        end
        
        plot_point_overlay(img, gmm, X);
        drawnow
        
        [gmm.n gmm.mu gmm.lam]

    end
  
end


