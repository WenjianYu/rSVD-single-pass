% Generate large matrix data using DCT.
m=2E5;                                 % size of matrix
n=2E5;
fid = fopen('exp1_2E5_2E5.dat','w');   % disk file name

% % % % % type 1% % % % %
S = zeros(n,1);
for i=1:20
    S(i)=10^(-4*(i-1)/19);
end

for i=21:n
   S(i)=(10^(-4))/(i-20)^(1/10);
end
% % % % % type 2% % % % %
%for i = 1:n
 %   S(i) = i^(-2);
%end

% % % % % type 3% % % % %
%for i = 1:n
 %   S(i) = i^(-3);
%end

s = spdiags(S,0,n,m);

for i=1:m
    E=zeros(m,1);
    E(i)=1;
    e=idct(E);
    f=s*e;
    F=idct(f);
    fwrite(fid, F,'float');
end
fclose(fid);
