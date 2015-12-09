function [US, TS, isC] = sortSchur(A, k)
    
    [U, T] = schur(A, 'real');
    es = ordeig(T);
    
    [~, ix] = sort(abs(es), 'descend');
    select = zeros(length(es), 1);
    select(ix(1:k)) = true;
    [US,TS] = ordschur(U, T, select);
    
    delta = (TS(k, k) - TS(k+1, k+1))^2 + 4 * TS(k+1, k) * TS(k, k+ ...
                                                      1);
    if delta < 0
        isC = 1;
    else
        isC = 0;
    end
    
end