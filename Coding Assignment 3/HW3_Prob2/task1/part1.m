
clear
close all

observed = readtable('Task1_obs100cycles_det.csv');
obs(:,1) = table2array(observed(:,2));
obs(:,2) = table2array(observed(:,3));

% l1_0 = 3.15;
% l2_0 = -0.4;
% l3_0 = 8;

%Adjusted based on a previous 30 iteration run

l1_0 = 3.264433593750001;
l2_0 = -0.380156250000000;
l3_0 = 8.078124999999996;

alph = 0.0;

params0 = [l1_0,l2_0,l3_0];

options = optimset('MaxFunEvals',30,'MaxIter',30,'PlotFcns',@optimplotfval);

params = fminsearch(@(params)simulation_error(obs, params, alph),params0,options);

outputTemp = python('main.py',num2str(params(1)),num2str(params(2)),num2str(params(3))); 
output = sscanf(outputTemp,'%f'); 

%% Plotting 

figure;
hold on
plot(output(1:100),'linewidth',1.5);
plot(obs(:,1),'linewidth',1.5);
title(['H | Weight = ', num2str(alph)]);
legend('Calibrated Model','Observed Data');
xlabel('Observation');
ylabel('Number of Patients at NonICU Hospital')
hold off;
grid on

figure;
hold on
plot(output(101:200),'linewidth',1.5);
plot(obs(:,2),'linewidth',1.5);
title(['ICU | Weight = ', num2str(1-alph)]);
legend('Calibrated Model','Observed Data');
xlabel('Observation');
ylabel('Number of Patients at ICU Hospital')
hold off;
grid on

figure;
histogram(output(1:100))
hold on
histogram(obs(:,1))
ylabel('Number of Patients at non-ICU Hospital')
title(['H | Weight = ', num2str(alph)]);


figure;
histogram(output(101:200))
hold on
histogram(obs(:,2))
ylabel('Number of Patients at ICU Hospital')
title(['ICU | Weight = ', num2str(1-alph)]);
