function [A, d]= PCAtestmatrix(m, n, t)
% function [A, d]= PCAtestmatrix(m, n, t)
% generate type t test matrix for PCA computing
% A is of mxn, and t= 1, 2, 3, 4, 5.
% d is the singular values of matrix A.
L = randn(m, m);
[U, ~] = qr(L);
L = randn(n, n);
[V, ~] = qr(L);
p= min(m,n);
d= zeros(p,1);
if t==1,
    d(1:20)= 10.^(-4/19*(0:19)');
    d(21:p)= 1e-4./((1:(p-20))'.^0.1);
elseif t==2,
    d= (1:p)'.^(-2);
elseif t==3,
    d= (1:p)'.^(-3);
elseif t==4,
    d= exp(-(1:p)'/7);
elseif t==5,
    d= 10.^(-0.1*(1:p)');
elseif t==6,
    d(1:3)= 1;
    d(4:6)= 0.67;
    d(7:9)= 0.34;
    d(10:12)= 0.01;
    d(13:p)= 1e-2*(p-13:-1:0)'./(p-13);
else
    disp('not-supportted t value');
end

S= spdiags(d, 0, m, n);
A = U * S * V;
