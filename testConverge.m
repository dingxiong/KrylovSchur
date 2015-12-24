function flag = testConverge(H, k, i, tol)
% test the convergence of the ith eigenvalue in Krylov-Schur iteration
% H is the hessenberg matrix after truncation
% The test rule is 
%    | b_i | < max( ||H(1:k, 1:k)||_F * epsilon, tol * | \lambda_i | ) 
% Parameters:
%    i          the ith position of vector b = H(k+1, 1:k)
%    k          the size of truncated H is (k+1) x k
%    epsilon    machine precesion
%    tol        tolerance passed by user
%    ||*||_F    Frobeneous norm
%    \lambda_i  the ith eigenvalue of H(1:k, 1:k). 
%                  1x1 block => \lambda_i = H(i, i)
%                  2x2 block => \lambda_i = eig( H(i:i+1, i:i+1))
%                  Here it is not possible for H(i-1:i, i-1:i) because
%                  whenever we find index i converges and it
%                  corresponds to a pair of complex eigenvalues,
%                  then we also mark i+1 converges.
% Return:
%    flag       -1 : not converge 
%                1 : real eigenvalue converges
%                2 : complex eigenvalue pair converges
% 
% Note :
%   if i < k, we test whether the eigenvalue is real or complex
%   if i = k, we assume ith eigenvalue is real.
    
    epsilon = 2e-16;
    if i < k
        delta = (H(i, i) - H(i+1, i+1))^2 + 4 * H(i+1, i) * H(i, i+ ...
                                                          1);
    else 
        delta = 1;
    end
    
    if delta > 0                        % real case
        if abs(H(k+1, i)) < max(norm(H, 'fro') * epsilon, abs(H(i, i)) * tol )
            flag = 1;
        else 
            flag = -1;
        end
        
    else                                % complex case
        lambda = (H(i,i) + H(i+1, i+1) + 1i*sqrt(-delta)) / 2;
        if abs(H(k+1, i)) < max(norm(H, 'fro') * epsilon, abs(lambda) * tol )
            flag = 2;
        else
            flag = -1;
        end
    end

end