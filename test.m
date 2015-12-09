global A;

N = 10;
k = 4;
m = 7;

A = rand(N, N);
v1 = rand(N, 1);
[Q, H] = KrylovSchur(@Ax, v1, N, k, m);

[ev, es] = eig(A); 
disp(['error of eig(A) : ', num2str(norm(A * ev - ev * es))]);
es = diag(es);
[~, ix] = sort(abs(es), 'descend');
es = es(ix);

disp(H);
disp(es);
disp('error of Krylov : ');
disp(A * Q(:, 1:k) - Q(:, 1:k+1) * H(1:k+1, 1:k) );
