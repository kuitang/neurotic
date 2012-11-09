function [ h ] = plot_intensity_hists( mm, X )
% Plot intensity histogram for each class.
    h = mm.h_diagnostic;            
    
    for kk = 1:mm.K        
        k = mm.k_idx(kk);
        subplot(2,mm.K,k);
        idxs = mm.s_z == k;
        x = X(idxs,3);
        hist(x);        
        title([' N = ' num2str(sum(idxs)) ...
               ' m = ' num2str(mean(x), 2) ...
               ' sd = ' num2str(std(x, 1), 2)]);   
        
        % GMM only. TODO: Generalize the parametric plotter.
        if k > 1
            subplot(2,mm.K,mm.K + k);
            pred_sd = sqrt(mm.pred_cov(3,3,k));
            pdf = @(x) tpdf( (x - mm.mean(k,3)) / pred_sd, mm.pred_dof(k) );
            ezplot(pdf, 0, 1);
            title(['fitted m = ' num2str(mm.mean(k,3), 2) ...
                   ' sd = ' num2str(pred_sd), 2]);
        end                
        
    end
        
end

