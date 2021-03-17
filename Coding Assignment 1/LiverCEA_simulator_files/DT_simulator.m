% Outputs: only keep two strategies on frontier, only cost and accuracy
% write to csv file
% Disease condition
% Pc = probability of cirrhosis in a given population
% Test characteristics
% senF4 = sensitivity of the testing method Fib4
% senM = sensitivity of the testing method MRE
% senFS = sensitivity of the testing method Fibroscan
% speF4 = specificity of the testing method Fib4
% speM = specificity of the testing method MRE
% speFS = specificity of the testing method Fibroscan
% Test cost
% C_Fib4 = cost of the testing method Fib4
% C_MRE = cost of the testing method MRE
% C_FS = cost of the testing method Fibroscan
%c_hcc,c_ev,c_bbevl,c_resect,c_transplant,
%conditional probability for EV and HCC
%condition-1 correctly diagnosed;2-misdiagnosed
%misdiagnose for hcc?
%pev_c = conditional probability P(EV|C)
% pev_nc = P(EV|NC)
% phcc= P(HCC|C)
%phcc= P(HCC|NC)
% test\condition    C         NC
%             +     CD1     MissC3
%             -     UD4       CR2
% This function is probably inside a markov process model.
function [result,strat,prob,cost]=DT_simulator(seed,c_bb,...
    c_hcc, C_MRE, C_FS, C_LBX, C_Fib4, p_bb, Pc, p_bld_lb,senHF4,senLF4, senM, senFS, ... 
    speM, speFS,senL,speL, undet,fmre,ffs, c_egd, c_ebl,  c_trans, c_rsct,c_bld_lb, speHcc, senHcc,...
    hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo,Tmax,pop)
%initialize all end nodes performance measures.
%Probabilities

leav=96;
% pop=zeros(leav,1);
prob=zeros(leav,1);
cost=zeros(leav,10);%6-facility;7-trans;8-rsct
eff=zeros(leav,4);
lb_bld=zeros(leav,2);
% CLDQ_non=1;
% (1-undet)=0.1;
% (1-undet)=0.6;
%=======================1. Fib4 + MRE
cost(1:16,1)=C_MRE+C_Fib4;
cost(3:4,1)=C_LBX+C_MRE+C_Fib4;
cost(7:8,1)=C_LBX+C_MRE+C_Fib4;
cost(11:12,1)=C_LBX+C_MRE+C_Fib4;
cost(15:16,1)=C_LBX+C_MRE+C_Fib4;
cost(17:18,1)=C_Fib4;
%1 c,Hfib4,+mre 1
prob(1,1)=(1-undet)*Pc*senHF4*senM*(1-fmre);
%2 c,Hfib4,-mre 4
prob(2,1)=(1-undet)*Pc*senHF4*(1-senM)*(1-fmre);
%3 c, +Hfib4,fmre, +lb,1
prob(3,1)=(1-undet)*Pc*senHF4*fmre*senL;
%4 c, +Hfib4,fmre, -lb,4
prob(4,1)=(1-undet)*Pc*senHF4*fmre*(1-senL);

%5 nc,+Hfib4,+mre 3
prob(5,1)=(1-undet)*(1-Pc)*(1-senLF4)*(1-speM)*(1-fmre);
%6 nc,+Hfib4,-mre 2
prob(6,1)=(1-undet)*(1-Pc)*(1-senLF4)*speM*(1-fmre);
%7 nc,+Hfib4,fmre,+lb ,3
prob(7,1)=(1-undet)*(1-Pc)*(1-senLF4)*fmre*(1-speL);
%8 nc,+Hfib4,fmre,-lb,2
prob(8,1)=(1-undet)*(1-Pc)*(1-senLF4)*fmre*speL;

%negative Hfib4
% %1 c,-Hfib4,+mre 1
% prob(9,1)=(1-undet)*(1-senHF4)*Pc*senM*(1-fmre);
% %2 c,-Hfib4,-mre 4
% prob(10,1)=(1-undet)*(1-senHF4)*Pc*(1-senM)*(1-fmre);
% %3 c, -Hfib4,fmre, +lb,1
% prob(11,1)=(1-undet)*(1-senHF4)*Pc*fmre*senL;
% %4 c, -Hfib4,fmre, -lb,4
% prob(12,1)=(1-undet)*(1-senHF4)*Pc*fmre*(1-senL);

% %5 nc,-Hfib4,+mre 3
% prob(13,1)=(1-undet)*speHF4*(1-Pc)*(1-speM)*(1-fmre);
% %6 nc,-Hfib4,-mre 2
% prob(14,1)=(1-undet)*speHF4*(1-Pc)*speM*(1-fmre);
% %7 nc,-Hfib4,fmre,+lb ,3
% prob(15,1)=(1-undet)*speHF4*(1-Pc)*fmre*(1-speL);
% %8 nc,-Hfib4,fmre,-lb,2
% prob(16,1)=(1-undet)*speHF4*(1-Pc)*fmre*speL;
% 


%1 c,Hfib4,+mre 1
prob(9,1)=undet*Pc*senM*(1-fmre);
%2 c,Hfib4,-mre 4
prob(10,1)=undet*Pc*(1-senM)*(1-fmre);
%3 c, +Hfib4,fmre, +lb,1
prob(11,1)=undet*Pc*fmre*senL;
%4 c, +Hfib4,fmre, -lb,4
prob(12,1)=undet*Pc*fmre*(1-senL);

%5 nc,+Hfib4,+mre 3
prob(13,1)=undet*(1-Pc)*(1-speM)*(1-fmre);
%6 nc,+Hfib4,-mre 2
prob(14,1)=undet*(1-Pc)*speM*(1-fmre);
%7 nc,+Hfib4,fmre,+lb ,3
prob(15,1)=undet*(1-Pc)*fmre*(1-speL);
%8 nc,+Hfib4,fmre,-lb,2
prob(16,1)=undet*(1-Pc)*fmre*speL;

%9 c,+Lfib4 4
prob(17,1)=(1-undet)*(1-senHF4)*Pc;
%9 c,-Lfib4 1
% prob(18,1)=(1-undet)*speLF4*Pc;
%10 nc,+Lfib4 2
prob(18,1)=(1-undet)*senLF4*(1-Pc);
%10 nc,-Lfib4 3
% prob(20,1)=(1-undet)*(1-senLF4)*(1-Pc);

%LBX bleeding probability and cost 
lb_bld(7:8,1)=p_bld_lb;
lb_bld(7:8,2)=c_bld_lb;
lb_bld(3:4,1)=p_bld_lb;
lb_bld(3:4,2)=c_bld_lb;
lb_bld(11:12,1)=p_bld_lb;
lb_bld(11:12,2)=c_bld_lb;
lb_bld(15:16,1)=p_bld_lb;
lb_bld(15:16,2)=c_bld_lb;
%===========================2. Fib4+LBX 
%REMEMBER: combine LBX facility costs
cost(19:26,1)=C_LBX+C_Fib4;
cost(27:28,1)=C_Fib4;
%15 c, +hfib4,+LB 1
prob(19,1)=(1-undet)*senHF4*Pc*senL;
%16 c, +hfib4,-LB 4
prob(20,1)=(1-undet)*senHF4*Pc*(1-senL);
%17 nc, +hfib4, +LB 3
prob(21,1)=(1-undet)*(1-senLF4)*(1-Pc)*(1-speL);
%18 nc, +hfib4, -LB 2
prob(22,1)=(1-undet)*(1-senLF4)*(1-Pc)*speL;

%21 c,ufib4,+LB 1
prob(23,1)=undet*Pc*senL;
%22 c,ufib4,-LB 4
prob(24,1)=undet*Pc*(1-senL);
%23 nc,ufib4,+LB 3
prob(25,1)=undet*(1-Pc)*(1-speL);
%24 nc,ufib4,-LB 2
prob(26,1)=undet*(1-Pc)*(speL);


%19 c,+Lfib4 4
prob(27,1)=(1-undet)*(1-senHF4)*Pc;
%19 c,-Lfib4 1
% prob(34,1)=(1-undet)*speLF4*Pc;
%20 nc,+Lfib4 2
prob(28,1)=(1-undet)*senLF4*(1-Pc);
%20 nc,-Lfib4 3
% prob(36,1)=(1-undet)*(1-senLF4)*(1-Pc);


lb_bld(19:26,1)=p_bld_lb;
lb_bld(19:26,2)=c_bld_lb;


%=============================3. Fib4+FS 
cost(29:44,1)=C_FS+C_Fib4;

cost(31:32,1)=C_LBX+C_FS+C_Fib4;
cost(35:36,1)=C_LBX+C_FS+C_Fib4;
cost(39:40,1)=C_LBX+C_FS+C_Fib4;
cost(43:44,1)=C_LBX+C_FS+C_Fib4;
cost(45:46,1)=C_Fib4;
%c,+hfib4,+fs 1
prob(29,1)=(1-undet)*senHF4*Pc*senFS*(1-ffs);
%c,+hfib4,-fs 4
prob(30,1)=(1-undet)*senHF4*Pc*(1-senFS)*(1-ffs);
%c, +hfib4,ffs, +lb,1
prob(31,1)=(1-undet)*senHF4*Pc*ffs*senL;
%c, +hfib4,ffs, -lb,4
prob(32,1)=(1-undet)*senHF4*Pc*ffs*(1-senL);

%nc,+hfib4,+fs 3
prob(33,1)=(1-undet)*(1-senLF4)*(1-Pc)*(1-speFS)*(1-ffs);
%nc,+hfib4,-fs 2
prob(34,1)=(1-undet)*(1-senLF4)*(1-Pc)*speFS*(1-ffs);
%nc,+hfib4,ffs,+lb ,3
prob(35,1)=(1-undet)*(1-senLF4)*(1-Pc)*ffs*(1-speL);
%nc,+hfib4,ffs,-lb,2
prob(36,1)=(1-undet)*(1-senLF4)*(1-Pc)*ffs*speL;

%c,+hfib4,+fs 1
prob(37,1)=undet*Pc*senFS*(1-ffs);
%c,+hfib4,-fs 4
prob(38,1)=undet*Pc*(1-senFS)*(1-ffs);
%c, +hfib4,ffs, +lb,1
prob(39,1)=undet*Pc*ffs*senL;
%c, +hfib4,ffs, -lb,4
prob(40,1)=undet*Pc*ffs*(1-senL);

%nc,+hfib4,+fs 3
prob(41,1)=undet*(1-Pc)*(1-speFS)*(1-ffs);
%nc,+hfib4,-fs 2
prob(42,1)=undet*(1-Pc)*speFS*(1-ffs);
%nc,+hfib4,ffs,+lb ,3
prob(43,1)=undet*(1-Pc)*ffs*(1-speL);
%nc,+hfib4,ffs,-lb,2
prob(44,1)=undet*(1-Pc)*ffs*speL;

% c, ufib4, +LB 1
% prob(61,1)=Pc*undet*senL;
% c, ufib4, -LB 4
% prob(62,1)=Pc*undet*(1-senL);
% nc,ufib4, +LB 3
% prob(63,1)=(1-Pc)*undet*(1-speL);
% nc, ufib4, -LB 2
% prob(64,1)=(1-Pc)*undet*speL;

%c,+lfib4 4
prob(45,1)=(1-undet)*(1-senHF4)*Pc;
%c,-lfib4 1
% prob(58,1)=(1-undet)*speLF4*Pc;
%nc,+lfib4 2
prob(46,1)=(1-undet)*senLF4*(1-Pc);
%nc,-lfib4 3
% prob(60,1)=(1-undet)*(1-senLF4)*(1-Pc);


%LBX bleeding probability and cost 
lb_bld(31:32,1)=p_bld_lb;
lb_bld(31:32,2)=c_bld_lb;
lb_bld(35:36,1)=p_bld_lb;
lb_bld(35:36,2)=c_bld_lb;
lb_bld(39:40,1)=p_bld_lb;
lb_bld(39:40,2)=c_bld_lb;
lb_bld(43:44,1)=p_bld_lb;
lb_bld(43:44,2)=c_bld_lb;


%===========================4. Fib4
cost(51:54)=C_LBX+C_Fib4;

% c,+hfib4 1
prob(47,1)=(1-undet)*senHF4*Pc;
% c,-hfib4 4
% prob(66,1)=(1-undet)*Pc * (1 - senHF4);
% nc,+hfib4 3
prob(48,1)=(1-undet)*(1-senLF4)*(1-Pc);
% nc,-hfib4 2
% prob(68,1)=(1-undet)*(1-Pc) * speHF4;

% c, +lfib4 4
prob(49,1)=(1-undet)*(1-senHF4)*Pc;
% c, -lfib4 1
% prob(70,1)=(1-undet)*Pc * speLF4;
% nc, +lfib4 2
prob(50,1)=(1-undet)*senLF4*(1-Pc);
% nc, -lfib4 3
% prob(72,1)=(1-undet)*(1-Pc) * (1-senLF4);


%c,+LB 1
prob(51,1)=undet*Pc * senL;
%c,-LB 4
prob(52,1)=undet*Pc * (1-senL);
%nc,+LB 3
prob(53,1)=undet*(1-Pc) * (1-speL);
%nc,-LB 2
prob(54,1)=undet*(1-Pc) * (speL);
%sum(prob(31:38))=1
lb_bld(51:54,1)=p_bld_lb;
lb_bld(51:54,2)=c_bld_lb;
%===============================5. MRE+LB
cost(55:58,1)=C_MRE+C_LBX;
cost(61:64,1)=C_MRE+C_LBX;
cost(59:60,1)=C_MRE;

%c, +MRE,+LB 1
prob(55,1)=Pc * senM * senL*(1-fmre);
%c,+MRE,-LB 4
prob(56,1)=Pc * senM*(1-senL)*(1-fmre);
%nc,+MRE,+LB 3
prob(57,1)=(1-Pc) * (1-speM)*(1-speL)*(1-fmre);
%nc,+MRE,-LB 2
prob(58,1)=(1-Pc) * (1-speM)*speL*(1-fmre);

%c, -MRE 4
prob(59,1)=Pc * (1-senM)*(1-fmre);
%nc,-MRE 2
prob(60,1)=(1-Pc)* speM*(1-fmre);

%c,fMRE,+LB 1
prob(61,1)=Pc * senL* fmre;
%c,fMRE,-LB 4
prob(62,1)=Pc * (1-senL)* fmre;
%nc,fMRE,+LB 3
prob(63,1)=(1-Pc) * (1-speL)* fmre;
%nc,fMRE,-LB 2
prob(64,1)=(1-Pc) * speL* fmre;


%sum(prob(39:44))=1

lb_bld(55:58,1)=p_bld_lb;
lb_bld(55:58,2)=c_bld_lb;
lb_bld(61:64,1)=p_bld_lb;
lb_bld(61:64,2)=c_bld_lb;
%=================================6. MRE
cost(69:72,1)=C_MRE+C_LBX;
cost(65:68,1)=C_MRE;

%c,+MRE 1
prob(65,1)=Pc * senM * (1-fmre); 
%c, -MRE 4
prob(66,1)=Pc * (1-senM)* (1-fmre);
%nc,+MRE 3  
prob(67,1)=(1-Pc) * (1-speM) * (1-fmre);
%nc, -MRE 2
prob(68,1)=(1-Pc)* speM * (1-fmre);

%c,fMRE,LB+, 1  
prob(69,1)=Pc * senL * fmre;
%c,fMRE,LB-, 4  
prob(70,1)=Pc * (1-senL) * fmre;
%nc, fMRE, LB+ 3
prob(71,1)=(1-Pc)* (1-speL) * fmre;
%nc, fMRE, LB- 2
prob(72,1)=(1-Pc)* speL * fmre;

lb_bld(69:72,1)=p_bld_lb;
lb_bld(69:72,2)=c_bld_lb;
%==========================7.FS+LB 
cost(73:76,1)=C_FS+C_LBX;
cost(79:82,1)=C_FS+C_LBX;
cost(77:78,1)=C_FS;

%c, +fs,+LB 1
prob(73,1)=Pc * senFS * senL*(1-ffs);
%c,+fs,-LB 4
prob(74,1)=Pc * senFS*(1-senL)*(1-ffs);
%nc,+fs,+LB 3
prob(75,1)=(1-Pc) * (1-speFS)*(1-speL)*(1-ffs);
%nc,+fs,-LB 2
prob(76,1)=(1-Pc) * (1-speFS)*speL*(1-ffs);

%c, -fs 4
prob(77,1)=Pc * (1-senFS)*(1-ffs);
%nc,-fs 2
prob(78,1)=(1-Pc)* speFS*(1-ffs);

%c,ffs,+LB 1
prob(79,1)=Pc * senL* ffs;
%c,ffs,-LB 4
prob(80,1)=Pc * (1-senL)* ffs;
%nc,ffs,+LB 3
prob(81,1)=(1-Pc) * (1-speL)* ffs;
%nc,ffs,-LB 2
prob(82,1)=(1-Pc) * speL* ffs;


%sum(prob(39:44))=1

lb_bld(73:76,1)=p_bld_lb;
lb_bld(73:76,2)=c_bld_lb;
lb_bld(79:82,1)=p_bld_lb;
lb_bld(79:82,2)=c_bld_lb;

%========================8.FS 
cost(87:90,1)=C_FS+C_LBX;
cost(83:86,1)=C_FS;

%c,+fs 1
prob(83,1)=Pc * senFS * (1-ffs); 
%c, -fs 4
prob(84,1)=Pc * (1-senFS)* (1-ffs);
%nc,+fs 3  
prob(85,1)=(1-Pc) * (1-speFS) * (1-ffs);
%nc, -fs 2
prob(86,1)=(1-Pc)* speFS * (1-ffs);

%c,ffs,LB+, 1  
prob(87,1)=Pc * senL * ffs;
%c,ffs,LB-, 4  
prob(88,1)=Pc * (1-senL) * ffs;
%nc, ffs, LB+ 3
prob(89,1)=(1-Pc)* (1-speL) * ffs;
%nc, ffs, LB- 2
prob(90,1)=(1-Pc)* speL * ffs;

lb_bld(87:90,1)=p_bld_lb;
lb_bld(87:90,2)=c_bld_lb;


%========================9. LB 
cost(91:94,1)=C_LBX;
%59 c, +LB 2
prob(91,1)=Pc*senL;
%60 c,-LB 4
prob(92,1)=Pc*(1-senL);
%61 nc,+LB 3
prob(93,1)=(1-Pc)*(1-speL);
%62 nc,-LB 2
prob(94,1)=(1-Pc)*speL;

lb_bld(91:94,1)=p_bld_lb;
lb_bld(91:94,2)=c_bld_lb;
%==================10 no tests===============
prob(95,1)=Pc;
prob(96,1)=1-Pc;
%================effectiveness assign
%%=============1c d
eff(1,1)=1;
eff(3,1)=1;
eff(9,1)=1;
eff(11,1)=1;
eff(19,1)=1;
eff(23,1)=1;
eff(29,1)=1;
eff(31,1)=1;
eff(37,1)=1;
eff(39,1)=1;
eff(47,1)=1;
eff(51,1)=1;
eff(55,1)=1;
eff(61,1)=1;
eff(65,1)=1;
eff(69,1)=1;
eff(73,1)=1;
eff(79,1)=1;
eff(83,1)=1;
eff(87,1)=1;
eff(91,1)=1;
%%=============2 c r
eff(6,2)=1;
eff(8,2)=1;
eff(14,2)=1;
eff(16,2)=1;
eff(18,2)=1;
eff(22,2)=1;
eff(26,2)=1;
eff(28,2)=1;
eff(34,2)=1;
eff(36,2)=1;
eff(42,2)=1;
eff(44,2)=1;
eff(46,2)=1;
eff(50,2)=1;
eff(54,2)=1;
eff(58,2)=1;
eff(60,2)=1;
eff(64,2)=1;
eff(68,2)=1;
eff(72,2)=1;
eff(76,2)=1;
eff(78,2)=1;
eff(82,2)=1;
eff(86,2)=1;
eff(90,2)=1;
eff(94,2)=1;
eff(96,2)=1;
%%=============3 misdiagnosed
eff(5,3)=1;
eff(7,3)=1;
eff(13,3)=1;
eff(15,3)=1;
eff(21,3)=1;
eff(25,3)=1;
eff(33,3)=1;
eff(35,3)=1;
eff(41,3)=1;
eff(43,3)=1;
eff(48,3)=1;
eff(53,3)=1;
eff(57,3)=1;
eff(63,3)=1;
eff(67,3)=1;
eff(71,3)=1;
eff(75,3)=1;
eff(81,3)=1;
eff(85,3)=1;
eff(89,3)=1;
eff(93,3)=1;
%%=============4 undiagnosed
eff(:,4)=ones(leav,1)-(eff(:,1)+eff(:,2)+eff(:,3));

% eff(2,4)=1;
% eff(4,4)=1;
% eff(9,4)=1;
% eff(12,4)=1;
% eff(16,4)=1;
% eff(19,4)=1;
% eff(22,4)=1;
% eff(26,4)=1;
% eff(28,4)=1;
% eff(33,4)=1;
% eff(36,4)=1;
% eff(40,4)=1;
% eff(44,4)=1;
% eff(48,4)=1;
% eff(51,4)=1;
% eff(54,4)=1;
% eff(58,4)=1;
% eff(62,4)=1;
% eff(66,4)=1;
% eff(69,4)=1;
% eff(72,4)=1;
% eff(76,4)=1;
% eff(80,4)=1;
% eff(84,4)=1;
% eff(87,4)=1;
%============================Simulation evaluation========================
% initialize result tables
num_res=21;
num_strat=10;
result=zeros(leav,num_res);
strat=zeros(num_strat,num_res);
rate=zeros(leav,6);
% QALY_CLDQ=zeros(leav,1);
%eff_sim=zeros(44,1);


%Simulate cohorts both correctly diagnosed or misdiagnosed
strat_list_cd=find(eff(:,1)==1);
strat_list_md=find(eff(:,3)==1);
condition_cd=ones(length(strat_list_cd),1) ;
condition_md=2*condition_cd;
strat_list=[strat_list_cd; strat_list_md];
condition=[condition_cd; condition_md];
l_list=length(strat_list);

for ii = 1:l_list
    cost_total=0;
    cost_ev_total=0;
    cost_hcc_total=0;
    cost_bb_total=0;
    cost_ebl_total=0;
    cost_rsct_total=0;
    cost_trans_total=0;
    bld_count=0;
    trans_count=0;
    rsct_count=0;
    death_ev_count=0;
    death_hcc_count=0;
    death_cp=0;
for jj = 1:pop  
    [seed,cost_jj,cost_egdjj,cost_hccjj,cost_bbjj,cost_ebljj,...
        cost_transjj,cost_rsctjj,death_evjj,death_hccjj,death_cpjj,bldjj,...
        surgeryjj]=sim_v4(seed,condition(ii),...
    Tmax,c_bb, c_hcc, p_bb,c_egd, c_ebl,c_trans, c_rsct,C_MRE,...
    speHcc, senHcc,hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo);
    cost_ev_total=cost_ev_total+cost_egdjj;
    cost_hcc_total=cost_hcc_total+cost_hccjj;
    cost_bb_total=cost_bb_total+cost_bbjj;
    cost_ebl_total=cost_ebl_total+cost_ebljj;
    cost_total=cost_total+cost_jj;
    cost_rsct_total=cost_rsct_total+cost_rsctjj;
    cost_trans_total=cost_trans_total+cost_transjj;
    bld_count=bld_count+bldjj;
    trans_count=trans_count+surgeryjj(1);
    rsct_count=rsct_count+surgeryjj(2);
    death_ev_count=death_ev_count+death_evjj;
    death_hcc_count=death_hcc_count+death_hccjj;
    death_cp=death_cp+death_cpjj;
end
cost(strat_list(ii),8)=cost_trans_total/pop;
cost(strat_list(ii),7)=cost_rsct_total/pop;
cost(strat_list(ii),6)=cost_ebl_total/pop;
cost(strat_list(ii),5)=cost_bb_total/pop;
cost(strat_list(ii),4)=cost_hcc_total/pop;
cost(strat_list(ii),3)=cost_ev_total/pop;
cost(strat_list(ii),2)=cost_total/pop;
rate(strat_list(ii),1)=bld_count;
rate(strat_list(ii),2)=trans_count;
rate(strat_list(ii),3)=rsct_count;
rate(strat_list(ii),4)=death_ev_count;
rate(strat_list(ii),5)=death_hcc_count;
rate(strat_list(ii),6)=death_cp;
end
% =========================================================



%Simulate for undiagnosed cohorts
strat_list2=find(eff(:,4)==1);
%condition2=ones(length(strat_list2),1);
l_list2=length(strat_list2);
%record notesting events
for ii = 1:l_list2
   death_ev_c=0;
   death_hcc_c=0;
   death_cp=0;
   bld_c=0;
   cost_fn=0;
%    QALY_total=0;
%    CLDQ_total=0;
   for jj = 1:pop
       [seed,costjj,death_ev,death_hcc,death_cpjj]=sim_v4_2(seed,...
    1,Tmax,c_bb, c_hcc, p_bb,c_egd, c_ebl,c_trans, c_rsct,C_MRE,...
    speHcc, senHcc,hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo);
   cost_fn=cost_fn+costjj;
   death_ev_c=death_ev_c+death_ev;
   death_hcc_c=death_hcc_c+death_hcc;
   death_cp=death_cp+death_cpjj;
   end
cost(strat_list2(ii),9)=cost_fn/pop;
rate(strat_list2(ii),1)=bld_c;
rate(strat_list2(ii),4)=death_ev_c;
rate(strat_list2(ii),5)=death_hcc_c;
rate(strat_list2(ii),6)=death_cp;
end
% ========================================================

% ======================================================
% results collection
%==================================================
result(:,1)=prob(:,1).*(cost(:,2)+cost(:,1)+cost(:,9)+lb_bld(:,1).*lb_bld(:,2));%total cost
result(:,2)=prob(:,1).*cost(:,1);
result(:,3)=prob(:,1).*cost(:,3);
result(:,4)=prob(:,1).*cost(:,4);
result(:,5)=prob(:,1).*cost(:,5);
result(:,6)=prob(:,1).*cost(:,6);
result(:,7)=prob(:,1).*cost(:,7);
result(:,8)=prob(:,1).*cost(:,8);

result(:,9)=prob(:,1).*lb_bld(:,1).*lb_bld(:,2);

result(:,10)=prob(:,1).*rate(:,1);
result(:,11)=prob(:,1).*rate(:,2);
result(:,12)=prob(:,1).*rate(:,3);
result(:,13)=prob(:,1).*rate(:,4);%ev death
result(:,14)=prob(:,1).*rate(:,5);%hcc death
result(:,15)=prob(:,1).*rate(:,6);%cp death

result(:,16)=prob(:,1).*(rate(:,4)+rate(:,5)+rate(:,6));%total death
result(:,17)=prob(:,1).*eff(:,1);%correctly diagnosed
result(:,18)=prob(:,1).*eff(:,2);%correctly ruled out
result(:,19)=prob(:,1).*eff(:,3);%
result(:,20)=prob(:,1).*eff(:,4);%


% Do nothing +
strat(1,:)=sum(result(1:18,:));

strat(2,:)=sum(result(19:28,:));

strat(3,:)=sum(result(29:46,:));

strat(4,:)=sum(result(47:54,:));

strat(5,:)=sum(result(55:64,:));

strat(6,:)=sum(result(65:72,:));

strat(7,:)=sum(result(73:82,:));

strat(8,:)=sum(result(83:90,:));

strat(9,:)=sum(result(91:94,:));

strat(10,:)=sum(result(95:96,:));

end
