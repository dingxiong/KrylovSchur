function [e, v, isC, flag, nc, ni] = KrylovSchurEig(Ax, v1, n, k, m, maxIt, tol)
% calculate eigenvalue and eigenvector of matrix A.
% For the description of paramters, please see function KrylovSchur().

    [Q, H, isC, flag, nc, ni] = KrylovSchur(Ax, v1, n, k, m, maxIt, tol);
    [V, D] = eig(H(1:k+isC, 1:k+isC));
    
    e = diag(D);
    v = Q(:, 1:k+isC) * V;
end