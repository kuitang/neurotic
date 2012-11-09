function [ samples ] = run_mm(iter_f, X, img, mm, niters, burn_in, sample_freq, recovery_file )    
    
    [N, D] = size(X);
    
    % Initialize
    samples = [];
    loglike = zeros(niters, 1);
    nz_idxs = mm.s_z > 0;    
    mm.n = accumarray(mm.s_z(nz_idxs), ones(1,sum(nz_idxs)));
    
    for k = 2:mm.K
        mm = gmm_recompute_cluster(mm, X, k);
    end
    
    mm.empty_clusters = (length(mm.n) + 1):mm.Kmax;
    remain = mm.Kmax - length(mm.n);
    mm.n = [ mm.n ; zeros(remain, 1) ];            
    
    % Parameters
    mm.mean = zeros(mm.Kmax, D);
    mm.cov  = zeros(D, D, mm.Kmax);
    mm.dof  = zeros(mm.Kmax, 1);
    mm.scale = zeros(mm.Kmax, 1);
    mm.pred_dof = zeros(mm.Kmax, 1);
    mm.pred_cov = zeros(D, D, mm.Kmax);
    mm.x_like   = zeros(N, mm.Kmax);
    
    mm.new_like = zeros(N, 1);
    for n = 1:N
        mm.new_like(n,:) = gmm_one_sample_posterior_like(mm, X(n,:));
    end 
    
    % The parameter mm.background_pdf works on scalars or vectors. We will
    % define mm.background_like to work on X.
    mm.background_like = @(mm, X) 1/mm.N * mm.background_pdf(X(:,3));
    
    % Precompute the Cholesky decomposition of prior_cov
    [mm.prior_stdev, p] = cholcov(mm.prior_cov);
    assert(p == 0, 'mm.prior_cov was not PSD!');
    
    % Setup the plotter
    h = figure;
    set(h, 'Units', 'normalized', 'position', [0.1 0.1 0.9 0.4]);
    
    % Handle crashes
    if nargin >=8 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end                
    
    for i = first:niters
        tic
        z_old = mm.s_z;
        
        % Iterate
        mm = iter_f(mm, X);        
        
        % Timing and diagnostics
        change = mean(abs(z_old ~= mm.s_z));
        
        itertime = toc;
        rate     = N / itertime;        
        
        fprintf(1, 'Iter %d: time = %0.2f; speed = %.2f; K = %d, change = %0.3f; loglike = %06.2f\n', ...
                i, itertime, rate, mm.K, change, mm.loglike);
        figure(10);
        nnz = mm.n(mm.n > 0);
        bar(nnz(2:end));        
        
        % Save the log-likelihood
        loglike(i) = mm.loglike;                
        
        %if mod(i, sample_freq) == 0
        if true
            % Draw the figure
            figure(h);
            subplot(1,3,1);
            plot(loglike(loglike ~= 0));
            title('Log-likelihood vs iteration');
            xlabel('Iteration');
            ylabel('Log-likelihood');
            
            subplot(1,3,2);
            plot_point_overlay(img, mm, X);
            title(['K = ' num2str(mm.K) ' iter ' num2str(i) ' ll = ' num2str(mm.loglike)]);
            
            subplot(1,3,3);
            imshow(img);            
            
            % Draw diagnostic images
            if isfield(mm, 'plot_diagnostic')
                % TODO: Proper handles
                figure(mm.h_diagnostic);
                mm.plot_diagnostic(mm, X);
            end            
            
            drawnow;
            
            if i > burn_in
                samples = [samples mm];
                iter = i;

                %save(['iter' num2str(i) '.mat']);
            end
        end        

    end
    
end


