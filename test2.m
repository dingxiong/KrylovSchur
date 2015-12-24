global A;

N = 10;
k = 4;
m = 7;

A = rand(N, N);
v1 = rand(N, 1);
[Q, H, isC] = KrylovSchur(@Ax, v1, N, k, m, 50);
es1 = ordeig(H);

[ev, es] = eig(A); 
disp(['error of eig(A) : ', num2str(norm(A * ev - ev * es))]);
es = diag(es);
[~, ix] = sort(abs(es), 'descend');
es = es(ix);

disp(es1);
disp(es);

