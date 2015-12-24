function [Q, H, isC, flag] = KrylovSchur(Ax, v1, n, k, m, maxIt)

    Q = zeros(n, m+1); 
    H = zeros(m+1, m);
    Q(:,1) = v1 / norm(v1);
    p = 1;                           % converge test position
                                        
    % initialize stage
    % At this stage, we do not need to conduct Schur decomposition 
    % because in the iteration stage, we will expand it and
    % decompose it.
    [Q, H] = expandKrylov(Ax, Q, H, 0, k); 
    % [U, T] = schur(H(1:k, 1:k), 'real');
    % H(1:k, 1:k) = T;
    % Q(:, 1:k) = Q(:, 1:k) * U;
    % H(k+1, 1:k) = H(k+1, k) * U(end, :);  % last row is enough
    
    isC = 0;
    % iteration stage
    for i = 1:maxIt
        if p > k, break; end
        
        % expand stage 
        [Q, H] = expandKrylov(Ax, Q, H, k+isC, m);
        [U, T, isC] = sortSchur(H(p:m, p:m), k-p+1); 
        %disp('err'); disp(Ax(Q(:, 1:m)) - Q(:, 1:m+1) * H(1:m+1,1:m)); 
        H(p:m, p:m) = T;
        H(1:p-1, p:m) = H(1:p-1, p:m) * U;
        Q(:, p:m) = Q(:, p:m) * U;
        H(m+1, p:m) = H(m+1, m) * U(end, :);    
        %disp('err'); disp(Ax(Q(:, 1:m)) - Q(:, 1:m+1) * H(1:m+1, 1:m)); 
        %disp('Q'); disp(Q); disp('H'); disp(H);
        
        % truncate stage
        [Q, H] = truncateKrylov(Q, H, k+isC, m); 
        %disp(Ax(Q(:, 1:k+isC)) - Q(:, 1:k+1+isC) * H(1:k+1+isC, 1:k+isC));
        %disp('Q'); disp(Q); disp('H'); disp(H);
        
        % test convergence
        check = true; disp(p);
        while check && p <= k
            result = testConverge(H, k+isC, p, 2e-16);
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
end