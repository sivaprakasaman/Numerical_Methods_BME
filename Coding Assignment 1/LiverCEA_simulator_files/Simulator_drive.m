% simulation driver script
% Relationship between the m files, and how to use them
% ==================================================
% Input arguments can be varied
% ==================================================

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

% ==================================================
%Termination time
Tmax=60;
%Population size
pop =100000;
seed=47906;
%main simulator
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
data_cells = num2cell(table);
% Header
row_header(1:3,1)={'Fib-4+MRE', 'Fib-4+LB', 'Fib-4+VCTE'};
col_header = {'Total cost', 'Accuracy'};
output_matrix=[{'Strategy'} col_header; row_header data_cells];     %Join cell arrays
%xlswrite(strcat('LiverCEA_results_',datestr(datetime('today')),'.xls'),output_matrix);     %Write data and both headers
writecell(output_matrix,strcat('LiverCEA_results_',datestr(datetime('today'))),"FileType",'spreadsheet');
