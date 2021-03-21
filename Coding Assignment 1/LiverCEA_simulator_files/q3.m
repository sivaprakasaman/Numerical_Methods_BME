%Andrew Sivaprakasam
%Coding Assignment 1 
%Question 3:


%% Compute Baselines:
% ==================================================
% Test characteristics
% ==================================================
 senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;
% ==================================================
% Disease condition
% ==================================================
Pc =.02;%Prevalence of cirrhosis
p_ev=.347; %prevalence of EV for people who already have cirrhosis  
p_evs=.48;  %prevalence of small EV if having EV
p_evto=ratecon(0.023,12); %Prevalence of developing EV, format(prevalence, duration)
hcc_inc=0.024;% hcc incidence
mor_hcco =[0.492 36]; %mortality if there was cancer, format: [prevalence duration]


p_bld_lb=0.006;% Bleeding probability for doing LB
p_bb=0.5; %probability of taking Beta blocker
% format for the following parameters: [prevalence duration]
r_bld_suto=[0.07 24]; %Bleeding prevalence if small EV
r_bld_luto=[0.25 24]; %Bleeding prevalence if large EV
r_bld_lbbo=[0.2 24]; %Bleeding prevalence if treated with Beta blocker
r_bld_leblo=[0.14 24];%Bleeding prevalence if treated with EBL
mor_bld=[0.163 6]; %mortality prevalence if there was bleeding
mor_cp=[0.02, 12]; %  mortality if having compensated cirrhosis
mor_transo=[0.17 36];%mortality if there was transplantation
mor_nohcco=[0.72 36];%mortality if there was no cancer
p_hccpreo=ratecon(0.025,12); %Prevalence of developing HCC

% ==================================================
% Test costs - baseline values
% ==================================================
c_bb=33.46; %beta block
c_hcc=136.76+22.85;%cost of afp is 22.85
C_MRE =544.18; %Cost of MRE
C_FS = 150.34; %Cost of FS
C_LBX=95.23+80+1236.62; %Cost of LB
C_Fib4=0 ; %Cost of Fib4 test
c_egd=699.79+73.4; %Cost of using EGD
c_ebl=1334.83+73.4; %Cost of using EBL
c_trans=739100; %Cost of transplantation
c_rsct=3156.5+18500; %cost of resection
c_bld_lb=1600; %Cost of treating LB caused bleeding

%% Running Simulation

%Termination time
Tmax=60;
%Population size
pop =1000;
seed=47906;

%Baseline
 [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

%Results collection 
num_strat = 3;
num_results = 2;
table = zeros(num_strat,num_results);
table(:,1) = stratbase(1:3,1);%total cost
table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr

cost_base(:,1) = table(:,1);
acc_base(:,1) = table(:,2);


%% Get full factorial design:

minmax = ff2n(10)+1;

%% max and min values for 10 parameters:
clear senL speL senHF4 senLF4 senM speM fmre sensFS speFS ffs

senL = [0.89,1]; %1
speL = [0.92,1]; %2
senHF4 = [0.35,.41]; %3
senLF4 = [0.74,.85]; %4
senM = [0.6,.97]; %5
speM = [0.84,.93]; %6
fmre = [0.04,.06]; %7
senFS = [0.78,.95]; %8
speFS = [0.85,.89]; %9
ffs = [0.035,.5]; %10

%% Run Sampling:

L = length(minmax);

Fib4MRE_cost_out = zeros(1,L);
Fib4LB_cost_out = zeros(1,L);
Fib4VCTE_cost_out = zeros(1,L);


Fib4MRE_acc_out = zeros(1,L);
Fib4LB_acc_out = zeros(1,L);
Fib4VCTE_acc_out = zeros(1,L);

tic 
parfor i = 1:L
    
    disp(['Iteration: ',num2str(i)]);
    table = zeros(num_strat,num_results);

        [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4(minmax(i,3)),senLF4(minmax(i,4)), senM(minmax(i,5)), senFS(minmax(i,8)), ...
        speM(minmax(i,6)), speFS(minmax(i,9)), senL(minmax(i,1)),speL(minmax(i,2)), undet,fmre(minmax(i,7)),ffs(minmax(i,10)), c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_out(i) = table(1,1);
    Fib4LB_cost_out(i) = table(2,1);
    Fib4VCTE_cost_out(i) = table(3,1);
    
    Fib4MRE_acc_out(i) = table(1,2);
    Fib4LB_acc_out(i) = table(2,2);
    Fib4VCTE_acc_out(i) = table(3,2);
    
end

toc

%% Plot the Sampled Distributions:

Fib4MRE_cost_mean = mean(Fib4MRE_cost_out);
Fib4MRE_cost_var = var(Fib4MRE_cost_out);

Fib4LB_cost_mean = mean(Fib4LB_cost_out);
Fib4LB_cost_var = var(Fib4LB_cost_out);

Fib4VCTE_cost_mean = mean(Fib4VCTE_cost_out);
Fib4VCTE_cost_var = var(Fib4VCTE_cost_out);


Fib4MRE_acc_mean = mean(Fib4MRE_acc_out);
Fib4MRE_acc_var = var(Fib4MRE_acc_out);

Fib4LB_acc_mean = mean(Fib4LB_acc_out);
Fib4LB_acc_var = var(Fib4LB_acc_out);

Fib4VCTE_acc_mean = mean(Fib4VCTE_acc_out);
Fib4VCTE_acc_var = var(Fib4VCTE_acc_out);

%dists_cost = figure;
figure;
subplot(1,3,1);
histogram(Fib4MRE_cost_out);
title('Fib4 + MRE Cost');
text(0,0,{['Mean: ', num2str(Fib4MRE_cost_mean)],['Var: ',num2str(Fib4MRE_cost_var)]});

subplot(1,3,2);
histogram(Fib4LB_cost_out);
title('Fib4 + LB Cost');
text(100,0,{['Mean: ', num2str(Fib4LB_cost_mean)],['Var: ',num2str(Fib4LB_cost_var)]});

subplot(1,3,3);
histogram(Fib4VCTE_cost_out);
title('Fib4 + VCTE Cost');
sgtitle('Distribution of Cost with Full Factorial Testing 10 Params')
text(200,0,{['Mean: ', num2str(Fib4VCTE_cost_mean)],['Var: ',num2str(Fib4VCTE_cost_var)]});


figure;
subplot(1,3,1);
histogram(Fib4MRE_acc_out);
title('Fib4 + MRE Accuracy');
text(1,0,{['Mean: ', num2str(Fib4MRE_acc_mean)],['Var: ',num2str(Fib4MRE_acc_var)]});

subplot(1,3,2);
histogram(Fib4LB_acc_out);
title('Fib4 + LB Accuracy');
text(1.0,0,{['Mean: ', num2str(Fib4LB_acc_mean)],['Var: ',num2str(Fib4LB_acc_var)]});

subplot(1,3,3);
histogram(Fib4VCTE_acc_out);
title('Fib4 + VCTE Accuracy');
sgtitle('Distribution of Accuracy with Full Factorial Testing 10 Params')
text(.9,0,{['Mean: ', num2str(Fib4VCTE_acc_mean)],['Var: ',num2str(Fib4VCTE_acc_var)]});


% dists_acc = figure;
% hold on
% histogram(Fib4MRE_acc_out);
% histogram(Fib4LB_acc_out);
% histogram(Fib4VCTE_acc_out);
% hold off
% title('Distribution of Accuracy with Full Factorial Testing 10 Params')
% legend('Fib4MRE','Fib4LB','Fib4VCTE')

%% 3a Computations:

%Run uncertainty analysis???



%Fitting Distributions to data:
pd_4MRE_cost = fitdist(Fib4MRE_cost_out','Normal');
pd_4LB_cost = fitdist(Fib4LB_cost_out','Normal');
pd_4VCTE_cost = fitdist(Fib4VCTE_cost_out','Normal');

pd_4MRE_acc = fitdist(Fib4MRE_acc_out','Normal');
pd_4LB_acc = fitdist(Fib4LB_acc_out','Normal');
pd_4VCTE_acc = fitdist(Fib4VCTE_acc_out','Normal');

figure;
subplot(2,1,1);
hold on
plot(pdf(pd_4MRE_cost,1:2000),"LineWidth",1.5)
plot(pdf(pd_4LB_cost,1:2000),"LineWidth",1.5)
plot(pdf(pd_4VCTE_cost,1:2000),"LineWidth",1.5)
legend('Fib4MRE','Fib4LB','Fib4VCTE')
title('Distribution of Cost with Full Factorial Testing 10 Params')
hold on

subplot(2,1,2);
hold on
plot(0:0.001:1.1,pdf(pd_4MRE_acc,0:0.001:1.1),"LineWidth",1.5)
plot(0:0.001:1.1,pdf(pd_4LB_acc,0:0.001:1.1),"LineWidth",1.5)
plot(0:0.001:1.1,pdf(pd_4VCTE_acc,0:0.001:1.1),"LineWidth",1.5)
hold off 
title('Distribution of Accuracy with Full Factorial Testing 10 Params')

%% 3b: Counting percentage of samples exceeding cost differential:

diff = abs(Fib4MRE_cost_out-Fib4LB_cost_out);
percent_greater = sum(diff>300)*100/L;

%% 3c: 10-fold full factorial sens analysis on MRE (senM, speM):
%refresh defaults 
senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;


senM = [0.6,.97]; 
speM = [0.84,.93]; 
senM = linspace(senM(1),senM(2),10);
speM = linspace(speM(1),speM(2),10);

design = fullfact([10,10]);

parfor i = 1:length(design)
    
     table = zeros(num_strat,num_results);

     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM(design(i,1)), senFS, ...
    speM(design(i,2)), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);


    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_compare1(i) = table(1,1);
    Fib4LB_cost_compare1(i) = table(2,1);
    
    Fib4MRE_acc_compare1(i) = table(1,2);
    Fib4LB_acc_compare1(i) = table(2,2);
    
end

diff = abs(Fib4MRE_cost_compare1-Fib4LB_cost_compare1);
diff_greater = (diff>300);
L = length(diff);
percent_greater_MRE_compare = sum(diff_greater)*100/L

%better visualization here:
for j = 1:L
    matr_1(design(j,2),design(j,1)) = diff_greater(j);
%     sens(j) = senM(design(j,1));
%     spec(j) = speM(design(j,2));
end

sens = senM(design(:,1));
spec = speM(design(:,2));

figure;
hold on;
scatter3(sens,spec,diff.*(diff<300),'sq','LineWidth',1.5);
scatter3(sens,spec,diff.*(diff>300),'sq','LineWidth',1.5);
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
title('Sensitivity Analysis of MRE sens/spec');
legend('Cost Diff. <300','Cost Diff. >300')
xlim([min(sens)-0.02, max(sens)+0.02]);
ylim([min(spec)-0.02, max(spec)+0.02]);
hold off

figure;
scatter(1:L,diff)
hold on 
plot(ones(L,1)*300,"LineWidth",2);
title('Sensitivity Analysis of MRE sens/spec')
ylabel('Cost ($)')
xlabel('Iteration');


%% 3d: 10-fold full factorial sens analysis on LB (senL, speL):
senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;


senL = [0.89,1]; %1
speL = [0.92,1]; %2 
senL = linspace(senL(1),senL(2),10);
speL = linspace(speL(1),speL(2),10);

design = fullfact([10,10]);

parfor i = 1:length(design)
    
     table = zeros(num_strat,num_results);

     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL(design(i,1)),speL(design(i,2)), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);


    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_compare2(i) = table(1,1);
    Fib4LB_cost_compare2(i) = table(2,1);
    
    Fib4MRE_acc_compare2(i) = table(1,2);
    Fib4LB_acc_compare2(i) = table(2,2);
    
end

diff = abs(Fib4MRE_cost_compare2-Fib4LB_cost_compare2);
diff_greater = (diff>300);
L = length(diff);
percent_greater_LB_compare = sum(diff_greater)*100/L

for j = 1:L
    matr_2(design(j,2),design(j,1)) = diff_greater(j);
end

sens = senL(design(:,1));
spec = speL(design(:,2));

figure;
hold on;
scatter3(sens,spec,diff.*(diff<300),'sq','LineWidth',1.5);
scatter3(sens,spec,diff.*(diff>300),'sq','LineWidth',1.5);
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
title('Sensitivity Analysis of MRE sens/spec');
legend('Cost Diff. <300','Cost Diff. >300')
xlim([min(sens)-0.02, max(sens)+0.02]);
ylim([min(spec)-0.02, max(spec)+0.02]);
hold off

figure;
scatter(1:L,diff)
hold on 
plot(ones(L,1)*300,"LineWidth",2);
title('Sensitivity Analysis of LB sens/spec')
ylabel('Cost ($)')
xlabel('Iteration');

%% 3e Sample-based Sensitivity Analysis on MRE (senM, speM):

senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;

senM = [0.6,.97]; 
speM = [0.84,.93];

random_samples = rand(20,2);
lhs_samples = lhsdesign(20,2);
orth_samples = ortho_sample(20,2);

senM_r = random_samples(:,1)*(senM(2)-senM(1))+senM(1);
speM_r = random_samples(:,2)*(speM(2)-speM(1))+speM(1);

senM_lhs = lhs_samples(:,1)*(senM(2)-senM(1))+senM(1);
speM_lhs = lhs_samples(:,2)*(speM(2)-speM(1))+speM(1);

senM_orth = orth_samples(:,1)*(senM(2)-senM(1))+senM(1);
speM_orth = orth_samples(:,2)*(speM(2)-speM(1))+speM(1);

figure;
subplot(1,3,1);
scatter(senM_r,speM_r);
hold on
plot(senM,ones(2,1)*mean(speM),'LineWidth',1.5)
plot(ones(2,1)*mean(senM),speM,'LineWidth',1.5)
xlim(senM)
ylim(speM)
xlabel('Sensitivity')
ylabel('Specificity')
title('Random');
axis square
hold off

subplot(1,3,2);
scatter(senM_lhs,speM_lhs);
hold on
plot(senM,ones(2,1)*mean(speM),'LineWidth',1.5)
plot(ones(2,1)*mean(senM),speM,'LineWidth',1.5)
xlim(senM)
ylim(speM)
xlabel('Sensitivity')
ylabel('Specificity')
title('LHS');
axis square
hold off

subplot(1,3,3);
scatter(senM_orth,speM_orth);
hold on
plot(senM,ones(2,1)*mean(speM),'LineWidth',1.5)
plot(ones(2,1)*mean(senM),speM,'LineWidth',1.5)
xlim(senM)
ylim(speM)
xlabel('Sensitivity')
ylabel('Specificity')
title('Orthogonal')
axis square
hold off

sgtitle('Verification of Sampling Distributions');


parfor i = 1:length(random_samples)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_r(i), senFS, ...
    speM_r(i), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_r(i) = table(1,1);
    Fib4LB_cost_r(i) = table(2,1);
    
    Fib4MRE_acc_r(i) = table(1,2);
    Fib4LB_acc_r(i) = table(2,2);
    
     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_lhs(i), senFS, ...
    speM_lhs(i), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_lhs(i) = table(1,1);
    Fib4LB_cost_lhs(i) = table(2,1);
    
    Fib4MRE_acc_lhs(i) = table(1,2);
    Fib4LB_acc_lhs(i) = table(2,2);
    
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_orth(i), senFS, ...
    speM_orth(i), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_orth(i) = table(1,1);
    Fib4LB_cost_orth(i) = table(2,1);
    
    Fib4MRE_acc_orth(i) = table(1,2);
    Fib4LB_acc_orth(i) = table(2,2);
    
    
end


mre_rand_reg_c = regress(Fib4MRE_cost_r',horzcat(senM_r,speM_r,ones(20,1)));
spe_corr_r_mre_cost = corr(Fib4MRE_cost_r',speM_r)
sen_corr_r_mre_cost = corr(Fib4MRE_cost_r',senM_r)
mre_lhs_reg_c = regress(Fib4MRE_cost_lhs',horzcat(senM_lhs,speM_lhs,ones(20,1)));
spe_corr_lhs_mre_cost = corr(Fib4MRE_cost_lhs',speM_lhs)
sen_corr_lhs_mre_cost = corr(Fib4MRE_cost_lhs',senM_lhs)
mre_orth_reg_c = regress(Fib4MRE_cost_orth',horzcat(senM_orth,speM_orth,ones(20,1)));
spe_corr_orth_mre_cost = corr(Fib4MRE_cost_orth',speM_orth)
sen_corr_orth_mre_cost = corr(Fib4MRE_cost_orth',senM_orth)


mre_rand_reg_a = regress(Fib4MRE_acc_r',horzcat(senM_lhs,speM_r,ones(20,1)));
spe_corr_r_mre_acc = corr(Fib4MRE_acc_r',speM_r)
sen_corr_r_mre_acc = corr(Fib4MRE_acc_r',senM_r)
mre_lhs_reg_a = regress(Fib4MRE_acc_lhs',horzcat(senM_lhs,speM_lhs,ones(20,1)));
spe_corr_lhs_mre_acc = corr(Fib4MRE_acc_lhs',speM_lhs)
sen_corr_lhs_mre_acc = corr(Fib4MRE_acc_lhs',senM_lhs)
mre_orth_reg_a = regress(Fib4MRE_acc_orth',horzcat(senM_orth,speM_orth,ones(20,1)));
spe_corr_orth_mre_acc = corr(Fib4MRE_acc_orth',speM_orth)
sen_corr_orth_mre_acc = corr(Fib4MRE_acc_orth',senM_orth)


figure;
plot3(senM_r,speM_r,Fib4MRE_cost_r,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Cost');
text(0.2,0.9,0,{['Cost = ',num2str(mre_rand_reg_c(1)),'senM + ', num2str(mre_rand_reg_c(2)),'speM + ', num2str(mre_rand_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_r_mre_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_r_mre_cost)]},'Units','Normalized')
title('Random Sampling')
saveas(gcf,'mre1','epsc')

figure;
plot3(senM_r,speM_r,Fib4MRE_acc_r,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Accuracy');
text(0.2,0.9,0,{['Accuracy = ',num2str(mre_rand_reg_a(1)),'senM + ', num2str(mre_rand_reg_a(2)),'speM + ', num2str(mre_rand_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_r_mre_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_r_mre_acc)]},'Units','Normalized')
title('Random Sampling')
saveas(gcf,'mre4','epsc')

figure;
plot3(senM_lhs,speM_lhs,Fib4MRE_cost_lhs,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Cost');
text(0.2,0.9,0,{['Cost = ',num2str(mre_lhs_reg_c(1)),'senM + ', num2str(mre_lhs_reg_c(2)),'speM + ', num2str(mre_lhs_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_lhs_mre_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_lhs_mre_cost)]},'Units','Normalized')
title('LH Sampling')
saveas(gcf,'mre2','epsc')

figure;
plot3(senM_lhs,speM_lhs,Fib4MRE_acc_lhs,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Accuracy');
title('LH Sampling')
text(0.2,0.9,0,{['Accuracy = ',num2str(mre_lhs_reg_a(1)),'senM + ', num2str(mre_lhs_reg_a(2)),'speM + ', num2str(mre_lhs_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_lhs_mre_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_lhs_mre_acc)]},'Units','Normalized')
saveas(gcf,'mre5','epsc')

figure;
plot3(senM_orth,speM_orth,Fib4MRE_cost_orth,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Cost');
text(0.2,0.9,0,{['Cost = ',num2str(mre_orth_reg_c(1)),'senM + ', num2str(mre_orth_reg_c(2)),'speM + ', num2str(mre_orth_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_orth_mre_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_orth_mre_cost)]},'Units','Normalized')
title('Orth Sampling')
saveas(gcf,'mre3','epsc')

figure;
plot3(senM_orth,speM_orth,Fib4MRE_acc_orth,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Accuracy');
title('Orth Sampling')
text(0.2,0.9,0,{['Accuracy = ',num2str(mre_orth_reg_a(1)),'senM + ', num2str(mre_orth_reg_a(2)),'speM + ', num2str(mre_orth_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_orth_mre_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_orth_mre_acc)]},'Units','Normalized')
saveas(gcf,'mre6','epsc')


 
%% 3e Sample-based Sensitivity Analysis on LB (senL, speL):

senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;

senL = [0.89,1]; %1
speL = [0.92,1]; %2

random_samples = rand(20,2);
lhs_samples = lhsdesign(20,2);
orth_samples = ortho_sample(20,2);

senL_r = random_samples(:,1)*(senL(2)-senL(1))+senL(1);
speL_r = random_samples(:,2)*(speL(2)-speL(1))+speL(1);

senL_lhs = lhs_samples(:,1)*(senL(2)-senL(1))+senL(1);
speL_lhs = lhs_samples(:,2)*(speL(2)-speL(1))+speL(1);

senL_orth = orth_samples(:,1)*(senL(2)-senL(1))+senL(1);
speL_orth = orth_samples(:,2)*(speL(2)-speL(1))+speL(1);


parfor i = 1:length(random_samples)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL_r(i),speL_r(i), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_r(i) = table(1,1);
    Fib4LB_cost_r(i) = table(2,1);
    
    Fib4MRE_acc_r(i) = table(1,2);
    Fib4LB_acc_r(i) = table(2,2);
    
     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL_lhs(i),speL_lhs(i), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_lhs(i) = table(1,1);
    Fib4LB_cost_lhs(i) = table(2,1);
    
    Fib4MRE_acc_lhs(i) = table(1,2);
    Fib4LB_acc_lhs(i) = table(2,2);
    
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL_orth(i),speL_orth(i), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_orth(i) = table(1,1);
    Fib4LB_cost_orth(i) = table(2,1);
    
    Fib4MRE_acc_orth(i) = table(1,2);
    Fib4LB_acc_orth(i) = table(2,2);
    
    
end


lb_rand_reg_c = regress(Fib4LB_cost_r',horzcat(senL_lhs,speL_r,ones(20,1)));
spe_corr_r_lb_cost = corr(Fib4LB_cost_r',speL_r)
sen_corr_r_lb_cost = corr(Fib4LB_cost_r',senL_r)
lb_lhs_reg_c = regress(Fib4LB_cost_lhs',horzcat(senL_lhs,speL_lhs,ones(20,1)));
spe_corr_lhs_lb_cost = corr(Fib4LB_cost_lhs',speL_lhs)
sen_corr_lhs_lb_cost = corr(Fib4LB_cost_lhs',senL_lhs)
lb_orth_reg_c = regress(Fib4LB_cost_orth',horzcat(senL_orth,speL_orth,ones(20,1)));
spe_corr_orth_lb_cost = corr(Fib4LB_cost_orth',speL_orth)
sen_corr_orth_lb_cost = corr(Fib4LB_cost_orth',senL_orth)


lb_rand_reg_a = regress(Fib4LB_acc_r',horzcat(senL_lhs,speL_r,ones(20,1)));
spe_corr_r_lb_acc = corr(Fib4LB_acc_r',speL_r)
sen_corr_r_lb_acc = corr(Fib4LB_acc_r',senL_r)
lb_lhs_reg_a = regress(Fib4LB_acc_lhs',horzcat(senL_lhs,speL_lhs,ones(20,1)));
spe_corr_lhs_lb_acc = corr(Fib4LB_acc_lhs',speL_lhs)
sen_corr_lhs_lb_acc = corr(Fib4LB_acc_lhs',senL_lhs)
lb_orth_reg_a = regress(Fib4LB_acc_orth',horzcat(senL_orth,speL_orth,ones(20,1)));
spe_corr_orth_lb_acc = corr(Fib4LB_acc_orth',speL_orth)
sen_corr_orth_lb_acc = corr(Fib4LB_acc_orth',senL_orth)

figure;
plot3(senL_r,speL_r,Fib4LB_cost_r,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
text(0.2,0.9,0,{['Cost = ',num2str(lb_rand_reg_c(1)),'senM + ', num2str(lb_rand_reg_c(2)),'speM + ', num2str(lb_rand_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_r_lb_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_r_lb_cost)]},'Units','Normalized')
title('Random Sampling')
saveas(gcf,'lb1','epsc')

figure;
plot3(senL_r,speL_r,Fib4LB_acc_r,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Accuracy');
text(0.2,0.9,0,{['Accuracy = ',num2str(lb_rand_reg_a(1)),'senM + ', num2str(lb_rand_reg_a(2)),'speM + ', num2str(lb_rand_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_r_lb_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_r_lb_acc)]},'Units','Normalized')
title('Random Sampling')
saveas(gcf,'lb4','epsc')

figure;
plot3(senL_lhs,speL_lhs,Fib4LB_cost_lhs,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
text(0.2,0.9,0,{['Cost = ',num2str(lb_lhs_reg_c(1)),'senM + ', num2str(lb_lhs_reg_c(2)),'speM + ', num2str(lb_lhs_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_lhs_lb_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_lhs_lb_cost)]},'Units','Normalized')
title('LH Sampling')
saveas(gcf,'lb2','epsc')


figure;
plot3(senL_lhs,speL_lhs,Fib4LB_acc_lhs,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
text(0.2,0.9,0,{['Accuracy = ',num2str(lb_lhs_reg_a(1)),'senM + ', num2str(lb_lhs_reg_a(2)),'speM + ', num2str(lb_lhs_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_lhs_lb_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_lhs_lb_acc)]},'Units','Normalized')
title('LH Sampling')
saveas(gcf,'lb5','epsc')


figure;
plot3(senL_orth,speL_orth,Fib4LB_cost_orth,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
text(0.2,0.9,0,{['Cost = ',num2str(lb_orth_reg_c(1)),'senM + ', num2str(lb_orth_reg_c(2)),'speM + ', num2str(lb_orth_reg_c(3))],['Corr Sen Cost = ', num2str(sen_corr_orth_lb_cost)]...
    ,['Corr Spe cost = ',num2str(spe_corr_orth_lb_cost)]},'Units','Normalized')
title('Orth Sampling')
saveas(gcf,'lb3','epsc')

figure;
plot3(senL_orth,speL_orth,Fib4LB_acc_orth,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
text(0.2,0.9,0,{['Accuracy = ',num2str(lb_orth_reg_a(1)),'senM + ', num2str(lb_orth_reg_a(2)),'speM + ', num2str(lb_orth_reg_a(3))],['Corr Sen Accuracy = ', num2str(sen_corr_orth_lb_acc)]...
    ,['Corr Spe Accuracy = ',num2str(spe_corr_orth_lb_acc)]},'Units','Normalized')
title('Orth Sampling')
saveas(gcf,'lb6','epsc')



%% 3H Sensitivity Analysis Using Samples with Induced Correlation: | MRE

% Rember to run this thrice!

correlation = 0.9;

senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;


z = lhsgeneral({makedist('Normal',mean(senM_lhs),std(senM_lhs)),makedist('Normal',mean(speM_lhs),std(speM_lhs))}, [1,correlation;correlation,1], 20);
senM_corr = z(:,1);
speM_corr = z(:,2);


parfor i = 1:length(random_samples)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_corr(i), senFS, ...
    speM_corr(i), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_corr(i) = table(1,1);
    Fib4LB_cost_corr(i) = table(2,1);
    
    Fib4MRE_acc_corr(i) = table(1,2);
    Fib4LB_acc_corr(i) = table(2,2);
       
end

figure;
plot3(senM_corr,speM_corr,Fib4MRE_cost_corr,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Cost');
title(['Correlated LHS Sampling | Correlation = ', num2str(correlation)])

figure;
plot3(senM_corr,speM_corr,Fib4MRE_acc_corr,'o','MarkerSize',3)
xlabel('MRE Sensitivity');
ylabel('MRE Specificity');
zlabel('MRE Accuracy');
title(['Correlated LHS Sampling | Correlation = ', num2str(correlation)])

spe_corr_corr_mre_cost = corr(Fib4MRE_cost_corr',speM_corr)
sen_corr_corr_mre_cost = corr(Fib4MRE_cost_corr',senM_corr)

spe_corr_corr_mre_acc = corr(Fib4MRE_acc_corr',speM_corr)
sen_corr_corr_mre_acc = corr(Fib4MRE_acc_corr',senM_corr)

%% 3H Sensitivity Analysis Using Samples with Induced Correlation: | LB
senHF4 = .38;%38
 senLF4 = .84;%84
 senM = .80;%80
 senFS = .80;%80
 speM = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;


z = lhsgeneral({makedist('Normal',mean(senL_lhs),std(senL_lhs)),makedist('Normal',mean(speL_lhs),std(speL_lhs))}, [1,correlation;correlation,1], 20);
senL_corr = z(:,1);
speL_corr = z(:,2);

parfor i = 1:length(random_samples)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
    speM, speFS,senL_corr(i),speL_corr(i), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_corr(i) = table(1,1);
    Fib4LB_cost_corr(i) = table(2,1);
    
    Fib4MRE_acc_corr(i) = table(1,2);
    Fib4LB_acc_corr(i) = table(2,2);
       
end

figure;
plot3(senL_corr,speL_corr,Fib4LB_cost_corr,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Cost');
title(['Correlated LHS Sampling | Correlation = ', num2str(correlation)])

figure;
plot3(senL_corr,speL_corr,Fib4LB_acc_corr,'o','MarkerSize',3)
xlabel('LB Sensitivity');
ylabel('LB Specificity');
zlabel('LB Accuracy');
title(['Correlated LHS Sampling | Correlation = ', num2str(correlation)])

spe_corr_corr_lb_cost = corr(Fib4LB_cost_corr',speL_corr)
sen_corr_corr_lb_cost = corr(Fib4LB_cost_corr',senL_corr)

spe_corr_corr_lb_acc = corr(Fib4LB_acc_corr',speL_corr)
sen_corr_corr_lb_acc = corr(Fib4LB_acc_corr',senL_corr)

%% 3I Importance Sampling of MRE: 

%Code adapted from http://www.inferencelab.com/importance-sampling-matlab-demo/
%Assuming the distribution of sensitivity and specificity in practice is
%actually non-uniform, rather it is biased towards the upper end,
%re-calculate mean cost/accuracy of Fib4 + MRE patients

% true probability distribution
true_func = @(x) betapdf(x,5,2);

% Do importance sampling
N = 20;
% uniform proposal distribution
x_samples = rand(N,1);
proposal = 1/N;
% evaluate for each sample
target = true_func(x_samples);
% calculate importance weight
w = target ./ proposal;
w = w ./ sum(w);
% resample, with replacement, according to importance weight
samples = randsample(x_samples,N,true,w);

senHF4 = .38;%38
 senLF4 = .84;%84
 senM_O = .80;%80
 senFS = .80;%80
 speM_O = .86;%86
 speFS = .81;%81
undet=0.32; %undetermined FS test results
 speHcc=.83;
 senHcc=.9;
 senL=.93;%93
 speL=.95;%95
 fmre=0.04;
 ffs=0.071;

senM = [0.6,.97]; 
speM = [0.84,.93];

senM_beta = samples(:,1)*(senM(2)-senM(1))+senM(1);
speM_beta = samples(:,1)*(speM(2)-speM(1))+speM(1);

% plot

figure;
subplot(2,2,1)
hist(senM_lhs,10)
title('LHS Sampling - Sensitivity')
axis square

subplot(2,2,2)
hist(speM_lhs,10)
title('LHS Sampling - Specificity')
axis square

subplot(2,2,3)
hist(senM_beta,10)
title('Beta Distribution - Sensitivity')
axis square

subplot(2,2,4)
hist(speM_beta,10)
title('Beta Distribution - Specificity')
axis square

%Varying Sensitivity 

parfor i = 1:length(senM_beta)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_beta(i), senFS, ...
    speM_O, speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_beta_sen(i) = table(1,1);    
    Fib4MRE_acc_beta_sen(i) = table(1,2);
       
end

%Varying Specificity
parfor i = 1:length(senM_beta)

     
     table = zeros(num_strat,num_results);
     
     [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM_O, senFS, ...
    speM_beta(i), speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop);

    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Fib4MRE_cost_beta_spe(i) = table(1,1);    
    Fib4MRE_acc_beta_spe(i) = table(1,2);
       
end

mean_cost_beta_sen = mean(Fib4MRE_cost_beta_sen);
mean_cost_beta_spe = mean(Fib4MRE_cost_beta_spe);

mean_acc_beta_sen = mean(Fib4MRE_acc_beta_sen);
mean_acc_beta_spe = mean(Fib4MRE_acc_beta_spe);



