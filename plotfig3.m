function plotfig3(A1, s1, A2, s2, k)
% function plotfig3(A1, s1, A2, s2, k)
% oversampling parameter is 10.
    subplot(1,2,1);
    [~, d1, ~]= rSVD_exSP(A1, k);
    d2 = rSVDbasic(A1, k);
    [~, d3, ~] = rSVDsp(A1, k);
    semilogy(1:k, d1, 1:k, d3(1:k), 1:k, d2, '--', 1:k, s1(1:k), '.', 'LineWidth', 1.5);
    legend('existing single-pass','our algorithm', 'basic randomized','SVD');
    ylabel('\sigma_{i}');
    
% below for siam resubmission!
%     semilogy(1:k, d1, 'x', 1:k, d3(1:k), 'ko',...
%             1:k, d2, '--', 1:k, s1(1:k), 'b-',  'LineWidth', 1, 'MarkerSize', 4);
%     legend('single-pass [2]','randQB\_FP', 'randQB','SVD');
%     (d1-s1(1:k)')./(d3(1:k)-s1(1:k)')
%     ylabel('\sigma_{j}');
    
    axis([1, 50, 1e-4, 1]);
    subplot(1,2,2);
    [~, d1, ~]= rSVD_exSP(A2, k);
    d2 = rSVDbasic(A2, k);
    [~, d3, ~] = rSVDsp(A2, k);
    semilogy(1:k, d1, 1:k, d3(1:k), 1:k, d2, '--', 1:k, s2(1:k), '.', 'LineWidth', 1.5);
    legend('existing single-pass','our algorithm', 'basic randomized','SVD');
    ylabel('\sigma_{i}');
    % below for siam resubmission!
%     semilogy(1:k, d1, 'x', 1:k, d3(1:k), 'ko',...
%             1:k, d2, '--', 1:k, s2(1:k), 'b-',  'LineWidth', 1, 'MarkerSize', 4);
%     legend('single-pass [2]','randQB\_FP', 'randQB','SVD');
%     (d1-s2(1:k)')./(d3(1:k)-s2(1:k)')
%     ylabel('\sigma_{j}');
    
    axis([1, 50, 1e-4, 1]);
end