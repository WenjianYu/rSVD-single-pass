function acc_valid(s2, smy, s1, d1, k)
% function acc_valid(s2, smy, s1, d1, k)
% Accuracy validation of computed k singular values and draw plots.
% need s2 (existing single-pass), smy (our single-pass), s1(basic randQB)
% d1 is the standard values of the singular values.
if nargin<5,
    k=50;
end
subplot(1,2,1);
semilogy(1:k, s2(1:k), 1:k, smy, 1:k, s1(1:k), '--', 1:k, d1(1:k), '.', 'LineWidth', 1.5);
axis([1, 50, 1e-5, 1]);
%axis([1, 50, 1e-4, 1]);
legend('existing single-pass','our algorithm', 'basic randomized','SVD');
ylabel('\sigma_{ii}');
subplot(1,2,2);
plot(1:k, s2(1:k)-d1(1:k), 1:k, smy-d1(1:k), 1:k, s1(1:k)-d1(1:k), '--', 'LineWidth', 1.5);
axis([1, 50, -2e-3, 14e-3]);
% axis([1, 50, -0.1, 0.75]);
ylabel('Error');
legend('existing single-pass','our algorithm', 'basic randomized');