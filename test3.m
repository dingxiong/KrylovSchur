% test eigenvalues
global A;

N = 100;
k = 5;
m = 10;

A = rand(N, N);

[ev, es] = eig(A); 
disp(['error of eig(A) : ', num2str(norm(A * ev - ev * es))]);
es = diag(es);
[~, ix] = sort(abs(es), 'descend');
es = es(ix);
disp('eigenvalues by matlab command eig() '); disp(es(1:k));

for i = 1:3
    fprintf(1, ' =======  trial i = %d  ========== \n', i);
    
    v1 = rand(N, 1);
    [e, v, isC, flag, nc, ni] = KrylovSchurEig(@Ax, v1, N, k, m, 200, 2e-16);
    disp(['error of KrylovSchurEig : ', num2str(norm(A * v - v * diag(e)))]);
    [~, ix] = sort(abs(e), 'descend');
    e = e(ix);
    
    fprintf(1, 'number of converged eigenvalues: %d \n', nc);
    fprintf(1, 'number of iteration used : %d \n', ni);
    disp('eigenvalues by KrylovSchurEig '); disp(e);
    
end
