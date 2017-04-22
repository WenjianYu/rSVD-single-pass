function [U, s, V]= rSVD_exSP(A, k, p, opt)
% function [U, s, V]= rSVD_exSP(A, k)
% The exisiting single-pass algorithm for randomized SVD.
if nargin <3
    p=10;
end
if nargin <4
    opt=1;
end
% p is for oversampling
% opt= 0 or 1.
[m, n]= size(A);
l= k+p;
Omg= randn(n, l);
Omgt= randn(m, l);
Y= A*Omg;
Yt= A'*Omgt;
[Q, R]= qr(Y, 0);
[Qt, Rt]= qr(Yt, 0);
if opt==0,
    B= (Omgt'*Q)\(Yt'*Qt);
    [Ut, S, Vt]= svd(B, 'econ');
    Ut= Q*Ut;
    U= Ut(:, 1:k);
    Vt= Qt*Vt;
    V= Vt(:, 1:k);
    s= diag(S);
    s= s(1:k);
elseif opt==1,
    B= (Omg'*Qt)\(Y'*Q);   % A' ~ Qt*B*Q' 
    [Ut, S, Vt]= svd(B, 'econ');
    Vt= Q*Vt;
    U= Vt(:, 1:k);
    Ut= Qt*Ut;
    V= Ut(:, 1:k);
    s= diag(S);
    s= s(1:k);
elseif opt==2,
    Q= Q(:, 1:k);
    Qt= Qt(:, 1:k);
    B= (Omgt'*Q)\(Rt(1:k, :)');
    [Ut, S, Vt]= svd(B, 'econ');
    U= Q*Ut;
    V= Qt*Vt;
    s= diag(S);  
else % opt=3
    Q= Q(:, 1:k);
    Qt= Qt(:, 1:k);
    B= (Omg'*Qt)\(R(1:k, :)');   % A' ~ Qt*B*Q' 
    [Ut, S, Vt]= svd(B, 'econ');
    U= Q*Vt;
    V= Qt*Ut;
    s= diag(S);
end

end