function [ m ] = struct_mean( st )
% m = struct_mean(st)
%
% stupid; there's a better way
    
    [~, N] = size(st);
    m = st(1);
    
    fields = fieldnames(st);   
    for i = 1:numel(fields)
        f = fields{i};
        % Delete non-numeric fields
        if ~isnumeric(st(1).(f))
            fields(i) = []
            i = i - 1;
        end
    end
    
    for n = 2:N
        for i = 1:numel(fields)            
            f = fields{i};            
            m.(f) = m.(f) + st(n).(f);            
        end            
    end
    
    for i = 1:numel(fields)
        f = fields{i};
        if isnumeric(m.(f))
            m.(f) = m.(f) / N;
        end
    end

end

