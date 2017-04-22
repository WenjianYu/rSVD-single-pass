% show different singular value distribution
% run PCAtestmatrix.m to generate the type 1~5 matrices firstly.
% suppose d1, d2, ..., d5 stores their singular values.
l= 50;
alld= [d1, d2, d3, d4, d5];
subplot(1,2,1);
plot((1:l)'*ones(1, 5),  alld(1:l,:), 'LineWidth', 1.5);
axis([1, 50, 0, 1]);
legend('Type1', 'Type2','Type3', 'Type4','Type5');
xlabel('(a)'); ylabel('\sigma_{ii}');
subplot(1,2,2);
semilogy((1:l)'*ones(1, 5), alld(1:l,:), 'LineWidth', 2);
axis([1, 50, 1e-6, 1]);
legend('Type1', 'Type2','Type3', 'Type4','Type5');
xlabel('(b)'); ylabel('\sigma_{ii}');

