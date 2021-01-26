%Andrew Sivaprakasam
%Numerical Methods
%Warmup 1
close all
clear 

load('arrays_to_sort.mat');
trials = 10; 
dims = [10, 20, 30, 50, 100, 500, 1000];

merge_times = zeros(trials, length(dims));
bubble_times = zeros(trials, length(dims));
insertion_times = zeros(trials, length(dims));
selection_times = zeros(trials, length(dims));

for i = 1:length(dims)
    
    switch dims(i) 
        case 10
            array = ten_dim;
        case 20
            array = twenty_dim;
        case 30
            array = thirty_dim;
        case 50
            array = fifty_dim;
        case 100
            array = hund_dim;
        case 500 
            array = randi(100,500,10);
        case 1000
            array = randi(100,1000,10);
    end
    
    for j = 1:trials
        
        to_sort = array(:,j);
        
        tic
        bubblesort(to_sort);
        bubble_times(j,i) = toc;
        
        tic
        insertionsort(to_sort);
        insertion_times(j,i) = toc;
        
        tic
        selectionsort(to_sort);
        selection_times(j,i) = toc;
        
        tic
        mergesort(to_sort);
        merge_times(j,i) = toc;
                       
    end

end

%% Plotting

hold on 

means = mean(merge_times);
stds = std(merge_times);

errorbar(dims, means, stds,'LineWidth',1.5);

means = mean(bubble_times);
stds = std(bubble_times);

errorbar(dims, means, stds,'LineWidth',1.5);

means = mean(selection_times);
stds = std(selection_times);

errorbar(dims, means, stds,'LineWidth',1.5);

means = mean(insertion_times);
stds = std(insertion_times);

errorbar(dims, means, stds,'LineWidth',1.5);

set(gca,'XTick',dims)
set(gca,'XTickLabel',dims)
set(gca,'XScale','log')
hold off

legend('MergeSort', 'BubbleSort', 'SelectionSort', 'InsertionSort','Location','Northwest');
xlabel('Dimension')
ylabel('Time (s)')
set(gca,'FontSize',12)

