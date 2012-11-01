function [ m ] = positive_to_interval( m_star )
% The inverse of m_star = 1 / (1 - m)
    m = m_star ./ (1 + m_star);

end

