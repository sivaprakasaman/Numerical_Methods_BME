%Andrew Sivaprakasam
%Numerical Methods
%Warmup 1

n = 1e3:1e3:50e3;
trials = 10000;

mult_times = zeros(trials,length(n));

for i = 1:trials

    for j = 1:length(n)

        tic
        multiply2(n,n);
        mult_times(i,j) = toc; 
        
    end

end

mult_avg = mean(mult_times);
mult_std = std(mult_times);

errorbar(n,mult_avg,mult_std,'k','LineWidth', 1.5);
title('Amortized Analysis of multiply')
xlabel('N');
ylabel('Runtime (s)')

set(gca,'FontSize',13);

figure;

subplot(1,2,1)
plot(n,mult_avg,'k','LineWidth',1.5);
title('Runtime of Multiply')
xlabel('N');
ylabel('Average Runtime (s)')
set(gca, 'FontSize', 12)

subplot(1,2,2)
plot(n,(n.*min(mult_times)).^log2(3),'k','LineWidth',1.5)
title('O(n^{log_b(a)})')
xlabel('N');
ylabel('Average Runtime (s)')

set(gca, 'FontSize', 12)