function [Q, H] = expandKrylov(Ax, Q, H, sk, ek)
% expand Krylov subspace.
%  The function contruct the sk+1, sk+2, ..., ek_th column of Q.
%  A * Q(:, 1:ek) = Q(:, 1:ek+1) * H
% Parameters:
%   sk       start index
%   ek       end index
% Return:
%   Q        the expanded orthornormal matrix with dimension [n x ek+1]
%   H        dimension [ek+1 x ek], the upper [ek x ek] block is Hessenberg

    for k = sk+1 : ek
        v = Ax(Q(:, k));                % A * Q(:, k)
        w = Q(:, 1:k)' * v;
        v = v - Q(:, 1:k) * w;
        w2 = Q(:, 1:k)' * v;            % double normalize
        v = v - Q(:, 1:k) * w2;
        w = w + w2;
        nv = norm(v); 
        Q(:, k+1) = v / nv;             % v_{k+1}
        H(1:k+1, k) = [w; nv];
    end

end