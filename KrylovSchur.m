function [Q,H] = KrylovSchur(Ax, v1, n, k, m)

    Q = zeros(n, m+1); 
    H = zeros(m+1, m);
    Q(:,1) = v1 / norm(v1);

    % initialize stage
    [Q, H] = expandKrylov(Ax, Q, H, 0, k); 
    [U, T] = schur(H(1:k, 1:k), 'real');
    H(1:k, 1:k) = U' * H(1:k, 1:k) * U;
    Q(:, 1:k) = Q(:, 1:k) * U;
    H(k+1, 1:k) = H(k+1, k) * U(end, :);  % last row is enough
    
    isC = 0;
    % iteration stage
    for i = 1:20
        % expand stage 
        [Q, H] = expandKrylov(Ax, Q, H, k+isC, m);
        [U, T, isC] = sortSchur(H(1:m, 1:m), k); % disp(Ax(Q(:, 1:m)) - Q(:, 1:m+1) * H(1:m+1, 1:m));
        H(1:m, 1:m) = T;
        Q(:, 1:m) = Q(:, 1:m) * U;
        H(m+1, 1:m) = H(m+1, m) * U(m, :);    %disp(Ax(Q(:, 1:m)) - Q(:, 1:m+1) * H(1:m+1, 1:m)); 
        %disp('Q'); disp(Q); disp('H'); disp(H);
        
        % truncate stage
        [Q, H] = truncateKrylov(Q, H, k+isC, m); %disp(Ax(Q(:, 1:k+isC)) - Q(:, 1:k+1+isC) * H(1:k+1+isC, 1:k+isC));
    end
end