function [ m ] = struct_mean( st )
% m = struct_mean(st)
%
% stupid; there's a better way
    
    [~, N] = size(st);
    m = st(1);
    
    fields = fieldnames(st);
    for n = 2:N
        for i = 1:numel(fields)
            f = fields{i};
            m.(f) = m.(f) + st(1).(f);            
        end
    end
    
    for i = 1:numel(fields)
        f = fields{i};
        m.(f) = m.(f) / N;
    end

end

