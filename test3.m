% test eigenvalues
global A;

N = 40;
k = 5;
m = 10;

A = rand(N, N);
v1 = rand(N, 1);
[e, v, isC, flag, nc, ni] = KrylovSchurEig(@Ax, v1, N, k, m, 100, 2e-16);
disp(['error of KrylovSchurEig : ', num2str(norm(A * v - v * diag(e)))]);

[ev, es] = eig(A); 
disp(['error of eig(A) : ', num2str(norm(A * ev - ev * es))]);
es = diag(es);
[~, ix] = sort(abs(es), 'descend');
es = es(ix);

disp('number of converged eigenvalues: '); disp(nc);
disp('number of iteration used ');   disp(ni);
disp('eigenvalues by KrylovSchurEig '); disp(e);
disp('eigenvalues by matlab command eig() '); disp(es(1:k));

