function [ samples ] = run_mm(iter_f, X, img, mm, niters, burn_in, sample_freq, recovery_file )    
    
    [N, D] = size(X);
    
    % Initialize
    samples = [];
    loglike = zeros(niters, 1);
    mm.n = accumarray(mm.s_z, ones(1,N));    
    
    % Setup the plotter
    h = figure;
    set(h, 'Units', 'normalized', 'position', [0.25 0.4 0.6 0.4]);
    
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
        
        fprintf(1, 'Iter %d: time = %0.2f; speed = %.2f; change = %0.3f; loglike = %06.2f\n', ...
                i, itertime, rate, change, mm.loglike);     
        
        % Save the log-likelihood
        loglike(i) = mm.loglike;                
        
        if mod(i, sample_freq) == 0                      
            % Draw the figure
            figure(h);
            subplot(1,2,1);
            plot(loglike(loglike ~= 0));
            title('Log-likelihood vs iteration');
            xlabel('Iteration');
            ylabel('Log-likelihood');
            
            subplot(1,2,2);
            plot_point_overlay(img, mm, X);
            title(['K = ' num2str(mm.K) ' iter ' num2str(i) ' ll = ' num2str(mm.loglike)]);
            
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


