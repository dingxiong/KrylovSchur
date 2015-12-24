function [US, TS, isC] = sortSchur(A, k)
% perform Schur decompostion on A and put the needed k eigenvalues
% on the upper left block.
% 
% Paramters:
%   A        a square matrix
%   k        number of eigenvalues that will be reordered
% 
% Return:
%   US       othormormal matrix
%   TS       quasi upper triangular matrix
%   isC      isC = 0, the kth eigenvalue is real. 
%            isC = 1, the kth and (k+1)th eigenvalues are complex
% 
% Note:
%   This function assumes that eigenvalues are reordered by their 
%   magnitudes. You can change this functions accordingly to obtain
%   eigenvalues that you need.
    
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