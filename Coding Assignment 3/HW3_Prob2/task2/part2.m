%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 3

clear
close all

observed = readtable('Task2_obs1cycle_sto.csv.csv');
obs(:,1) = table2array(observed(:,1));
obs(:,2) = table2array(observed(:,2));

l1_0 = 3.15;
l2_0 = -0.4;
l3_0 = 8;

seed = table2array(readtable('seeds100.csv','ReadVariableNames',false));

output1 = zeros(30,2);

parfor i = 1:30
    outputTemp = python('main.py','3.15','-0.4','8',num2str(seed(i))); %e.g. l1:'3.15', l2:'-0.2', l3:'8' 
    output1(i,:) = transpose(sscanf(outputTemp,'%f')); %
end

% qqplot(output1(:,1))
% title('QQ Plot of H - Model Data')
% figure;
% qqplot(output1(:,2))
% title('QQ Plot of ICU - Model Data')
% 
% 
% qqplot(obs(:,1))
% title('QQ Plot of H - Observed Data')
% figure;
% qqplot(obs(:,2))
% title('QQ Plot of ICU - Observed Data')

%compare initial distributions between modeled and observed:
figure;
histogram(output1(:,1))
hold on
histogram(obs(:,1))
ylabel('Number of Beds at non-ICU Hospital')
title('H | Uncalibrated');
legend('Modeled','Observed')
xlim([130,190])

figure;
histogram(output1(:,2))
hold on
histogram(obs(:,2))
ylabel('Number of Beds at ICU Hospital')
title('ICU | Uncalibrated');
legend('Modeled','Observed')
xlim([60,110])

%% Calibration
alph = 0.5;

params0 = [l1_0,l2_0,l3_0];

options = optimset('MaxFunEvals',50,'MaxIter',50,'PlotFcns',@optimplotfval);

params = fminsearch(@(params)simulation_error(obs, params, seed, alph),params0,options);

output = zeros(30,2);

l1 = num2str(params(1));
l2 = num2str(params(2));
l3 = num2str(params(3));

parfor i = 1:30
    outputTemp = python('main.py',l1,l2,l3,num2str(seed(i))); %e.g. l1:'3.15', l2:'-0.2', l3:'8' 
    output(i,:) = transpose(sscanf(outputTemp,'%f')); %
end


%% Final plots 
figure;
histogram(output(:,1))
hold on
histogram(obs(:,1))
ylabel('Number of Beds at non-ICU Hospital')
title('H | Calibrated');
legend('Modeled','Observed')
xlim([130,190])


figure;
histogram(output(:,2))
hold on
histogram(obs(:,2))
ylabel('Number of Beds at ICU Hospital')
title('ICU | Calibrated');
legend('Modeled','Observed')
xlim([60,110])

