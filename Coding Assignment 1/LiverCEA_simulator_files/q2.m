%Andrew Sivaprakasam
%Coding Assignment 1 
%Question 2:

%% Default Parameters:
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

%Vary Pc
minmax = [0.0027,.02];
num_strat = 3;
num_results = 2;
table = zeros(num_strat,num_results);

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
    
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, minmax(i), p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    Pc_cost_out(:,i) = table(:,1);
    Pc_acc_out(:,i) = table(:,2);
end


%% HCC

table = zeros(num_strat,num_results);

minmax = [0.464,77.2];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,[minmax(i) mor_hcco(2)],mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    hcc_cost_out(:,i) = table(:,1);
    hcc_acc_out(:,i) = table(:,2);
end

%% Liver Biopsy Sensitivity

table = zeros(num_strat,num_results);

minmax = [0.89,1];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS, minmax(i),speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    senL_cost_out(:,i) = table(:,1);
    senL_acc_out(:,i) = table(:,2);
end

%% LB Spec
table = zeros(num_strat,num_results);

minmax = [0.92,1];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS, senL,minmax(i), undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    speL_cost_out(:,i) = table(:,1);
    speL_acc_out(:,i) = table(:,2);
end

%% senHF4

table = zeros(num_strat,num_results);

minmax = [0.35,.41];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,minmax(i),senLF4, senM, senFS, ...
        speM, speFS, senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    senHF4_cost_out(:,i) = table(:,1);
    senHF4_acc_out(:,i) = table(:,2);
end

%% senLF4
table = zeros(num_strat,num_results);
minmax = [0.74,.85];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,minmax(i), senM, senFS, ...
        speM, speFS, senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    senLF4_cost_out(:,i) = table(:,1);
    senLF4_acc_out(:,i) = table(:,2);
end

%% senM

table = zeros(num_strat,num_results);

minmax = [0.6,.97];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, minmax(i), senFS, ...
        speM, speFS, senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    senM_cost_out(:,i) = table(:,1);
    senM_acc_out(:,i) = table(:,2);
end

%% speM


table = zeros(num_strat,num_results);

minmax = [0.84,.93];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        minmax(i), speFS, senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    speM_cost_out(:,i) = table(:,1);
    speM_acc_out(:,i) = table(:,2);
end

%% fmre


table = zeros(num_strat,num_results);

minmax = [0.04,.06];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS, senL,speL, undet,minmax(i),ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    fmre_cost_out(:,i) = table(:,1);
    fmre_acc_out(:,i) = table(:,2);
end

%% senFs

table = zeros(num_strat,num_results);

minmax = [0.78,.95];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, minmax(i), ...
        speM, speFS, senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    senFs_cost_out(:,i) = table(:,1);
    senFs_acc_out(:,i) = table(:,2);
end

%% speFs

table = zeros(num_strat,num_results);

minmax = [0.85,.89];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, minmax(i), senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    speFs_cost_out(:,i) = table(:,1);
    speFs_acc_out(:,i) = table(:,2);
end

%% ffs

table = zeros(num_strat,num_results);

minmax = [0.035,.5];
%main simulator
for i = 1:length(minmax) 
    [result,stratbase,prob,cost]=DT_simulator(seed,c_bb,...
        c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ...
        speM, speFS, senL,speL, undet,fmre,minmax(i), c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
        hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
    mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
    p_hccpreo,Tmax,pop);
    table(:,1) = stratbase(1:3,1);%total cost
    table(:,2) = stratbase(1:3,17)+stratbase(1:3,18);%cd+cr
    
    %order 'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'
    ffs_cost_out(:,i) = table(:,1);
    ffs_acc_out(:,i) = table(:,2);
end

%% Plotting:

%'Fib-4+MRE'

%Cost
r = 1;
c = 1;
ticklabels = {"Pc", "mor_hcco", "senL", "speL", "senHF4", "senLF4","senM","speM","fmre","senFS","speFs","ffs"};
cost_min = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_min = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];

c = 2;
cost_max = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_max = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];

bar_length_cost = abs(cost_min-cost_base(r))+abs(cost_max-cost_base(r));
bar_length_acc = abs(acc_min-acc_base(r))+abs(acc_max-acc_base(r));

[~,I_c] = sort(bar_length_cost);
[~,I_a] = sort(bar_length_acc);
%min
figure;
subplot(1,2,1);
hold on
barh(cost_min(I_c)','BaseValue', cost_base(r))
barh(cost_max(I_c)','BaseValue', cost_base(r))
yticks([1:12])
yticklabels(ticklabels(I_c));

%max
subplot(1,2,2);
hold on
barh(acc_min(I_a)','BaseValue', acc_base(r))
barh(acc_max(I_a)','BaseValue', acc_base(r))
legend("Min","Max",'Location','southeast');
sgtitle("Fib-4+MRE");
yticks([1:12])
yticklabels(ticklabels(I_a));

%'Fib-4+LB'


%Cost
r = 2;
c = 1;
ticklabels = {"Pc", "mor_hcco", "senL", "speL", "senHF4", "senLF4","senM","speM","fmre","senFS","speFs","ffs"};
cost_min = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_min = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];

c = 2;
cost_max = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_max = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];


bar_length_cost = abs(cost_min-cost_base(r))+abs(cost_max-cost_base(r));
bar_length_acc = abs(acc_min-acc_base(r))+abs(acc_max-acc_base(r));

[~,I_c] = sort(bar_length_cost);
[~,I_a] = sort(bar_length_acc);

%min
figure;
subplot(1,2,1);
hold on
barh(cost_min(I_c)','BaseValue', cost_base(r))
barh(cost_max(I_c)','BaseValue', cost_base(r))
yticks([1:12])
yticklabels(ticklabels(I_c));

%max
subplot(1,2,2);
hold on
barh(acc_min(I_a)','BaseValue', acc_base(r))
barh(acc_max(I_a)','BaseValue', acc_base(r))
legend("Min","Max",'Location','southeast');
sgtitle("Fib-4+LB");
yticks([1:12])
yticklabels(ticklabels(I_a));


%'Fib-4+VCTE'


%Cost
r = 3;
c = 1;
ticklabels = {"Pc", "mor_hcco", "senL", "speL", "senHF4", "senLF4","senM","speM","fmre","senFS","speFs","ffs"};
cost_min = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_min = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];

c = 2;
cost_max = [Pc_cost_out(r,c),hcc_cost_out(r,c),senL_cost_out(r,c),speL_cost_out(r,c),senHF4_cost_out(r,c),senLF4_cost_out(r,c),senM_cost_out(r,c),speM_cost_out(r,c),fmre_cost_out(r,c),senFs_cost_out(r,c),speFs_cost_out(r,c),ffs_cost_out(r,c)];
acc_max = [Pc_acc_out(r,c),hcc_acc_out(r,c),senL_acc_out(r,c),speL_acc_out(r,c),senHF4_acc_out(r,c),senLF4_acc_out(r,c),senM_acc_out(r,c),speM_acc_out(r,c),fmre_acc_out(r,c),senFs_acc_out(r,c),speFs_acc_out(r,c),ffs_acc_out(r,c)];


bar_length_cost = abs(cost_min-cost_base(r))+abs(cost_max-cost_base(r));
bar_length_acc = abs(acc_min-acc_base(r))+abs(acc_max-acc_base(r));

[~,I_c] = sort(bar_length_cost);
[~,I_a] = sort(bar_length_acc);


%min
figure;
subplot(1,2,1);
hold on
barh(cost_min(I_c)','BaseValue', cost_base(r))
barh(cost_max(I_c)','BaseValue', cost_base(r))
yticks([1:12])
yticklabels(ticklabels(I_c));

%max
subplot(1,2,2);
hold on
barh(acc_min(I_a)','BaseValue', acc_base(r))
barh(acc_max(I_a)','BaseValue', acc_base(r))
legend("Min","Max",'Location','southeast');
sgtitle("Fib-4+VCTE");
yticks([1:12])
yticklabels(ticklabels(I_a));

