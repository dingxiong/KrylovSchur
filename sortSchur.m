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
%   
%  Also, take care when to determine isC, which should be
%  determined by T not TS.
    
    [U, T] = schur(A, 'real');
    es = ordeig(T);     
    [~, ix] = sort(abs(es), 'descend');
    
    % judge the k-th and (k+1)-th eigenvalues are conjugate complex
    % pair or not.
    k1 = ix(k);
    k2 = ix(k+1);
    % here must use k2 not k1+1 because k1+1 maybe out of bounds
    delta = (T(k1, k1) - T(k2, k2))^2 + 4 * T(k2, k1) * T(k1, k2);
    if k2 - k1 == 1 && delta < 0
        isC = 1;
    else 
        isC = 0;
    end
    
    % the following method has numerical round-off issue
    %e1 = es(ix(k));
    %e2 = es(ix(k+1));
    %isC = ~isreal(e1) && isreal(e1+e2);
    
    select = zeros(length(es), 1);
    select(ix(1:k+isC)) = true;
    [US, TS] = ordschur(U, T, select); 
    
    %disp([(1:size(es, 1))', abs(es), ix]); 
    %disp([(1:k+isC)', abs(ordeig(TS(1:k+isC, 1:k+isC)))]);
    
end