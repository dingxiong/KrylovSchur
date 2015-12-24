function [Q, H] = truncateKrylov(Q, H, k, m)
% truncate Krylov subsapce
    
%   Q(:, k+1) = Q(:, m+1);
%   H(k+1, 1:k) = H(m+1, 1:k);
%   H(m+1, 1:m) = zeros(1, m);               % reset the written part to zeros
                                             
   Q = [Q(:, 1:k), Q(:, m+1)];
   H = [H(1:k, 1:k); H(m+1, 1:k)];  % disp('Q'); disp(Q); disp('H'); disp(H);
end                                     