% for drawing Fig. 5 of our IJCAI paper.
% U0 is the standard result; Umy is from our single-pass algorithm.
[U0sort, IX]= sort(U0(:,1));

maxerr= zeros(50,1);
coeff= zeros(50,1);
for i=1:50;
    maxerr(i)= norm(U0(:,i)+Umy(:,i), 'inf');
    corref= corrcoef(U0(:,i),Umy(:,i));
    coeff(i)= abs(corref(1,2));
end
subplot(1,2,1);
plot(1:3000, U0sort, 1:3000, -Umy(IX,1), '+');
xlabel('(a)');
legend('SVD', 'Our alg.');
subplot(1,2,2);
plot(1:10, coeff(1:10), 'o-');
% plot(1:50, maxerr);
axis([1,10, 0.99, 1]);
xlabel('(b)');
ylabel('Correlation coefficient');