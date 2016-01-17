function [Q, H, isC, flag, nc, ni] = KrylovSchur(Ax, v1, n, k, m, maxIt, tol)
% conduct Krylov-Schur decomposition to find invariant subspace of
% matrix A.
% A * Q(:, 1:k+isC) = Q(:, 1:k+isC) * H(1:k+isC, 1:k+isC)
% =======================================================
% Reference:
% A KRYLOV–SCHUR ALGORITHM FOR LARGE EIGENPROBLEMS
% G. W. STEWART
% SIAM J. M ATRIX A NAL. A PPL .
% Vol. 23, No. 3, pp. 601–614
% =======================================================
% Paramters:
%   Ax       function takes vector x and returns A*x
%   v1       initial vector
%   n        size of matrix A
%   k        number of eigenvalues needed
%   m        Krylov restart threshold (largest dimension of krylov subspace)
%   maxIt    maximum iteration number
%   tol      convergence tolerance
% 
% Return:
%   Q        orthonormal matrix with dimension [n x k+1] or [n x k+2]
%   H        `Hessenberg' matrix with dimension [k+1 x k] or [k+2 x k+1]
%   isC      isC = 0, the kth eigenvalue is real. 
%            isC = 1, the kth and (k+1)th eigenvalues are complex
%            conjugate pair, so Q and H has one more dimension
%  flag      flag = 0, converge.
%            flat = 1, not converge
%  nc        number of converged eigenvalues 
%  ni        number of iteration used (each whole expansion stage
%            counts 1)

    Q = zeros(n, m+1); 
    H = zeros(m+1, m);
    Q(:,1) = v1 / norm(v1);
    p = 1;                              % converge test position
    isC = 0;                            % complex bit
                                        
    % initialize stage
    % At this stage, we do not need to conduct Schur decomposition 
    % because in the iteration stage, we will expand it and
    % decompose it.
    [Q, H] = expandKrylov(Ax, Q, H, 0, k); 
    
   
    % iteration stage
    i = 0;
    while i < maxIt && p <= k
        i = i+1;
        fprintf(1, 'i = %d, isC = %d, p = %d, r = %g \n', ...
                i, isC, p, H(k+1+isC, p));
        
        % expand stage 
        [Q, H] = expandKrylov(Ax, Q, H, k+isC, m);
        [U, T, isC] = sortSchur(H(p:m, p:m), k-p+1); 
        H(p:m, p:m) = T;
        H(1:p-1, p:m) = H(1:p-1, p:m) * U;
        Q(:, p:m) = Q(:, p:m) * U;
        H(m+1, p:m) = H(m+1, m) * U(end, :);    
        %disp('err'); disp(Ax(Q(:, 1:m)) - Q(:, 1:m+1) * H(1:m+1, 1:m)); 
        %disp('Q'); disp(Q); disp('H'); disp(H);
        
        % truncate stage
        [Q, H] = truncateKrylov(Q, H, k+isC, m); 
        % disp(Ax(Q(:, 1:k+isC)) - Q(:, 1:k+1+isC) * H(1:k+1+isC, 1:k+isC));
        
        % test convergence
        check = true; 
        while check
            result = testConverge(H, k+isC, p, tol);
            if result == 1 || result == 2
                p = p + result;
                if p > k
                    check = false;
                end
            else
                check = false;
            end
        end
        
    end
    
    % return the convergence infomation
    ni = i;
    if p > k
        flag = 0;                       % converges
        nc = k + isC;
    else
        flag = 1;                       % not converge
        nc = p - 1;
    end
        
end