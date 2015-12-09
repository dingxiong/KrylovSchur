function [Q, H] = expandKrylov(Ax, Q, H, sk, ek)

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