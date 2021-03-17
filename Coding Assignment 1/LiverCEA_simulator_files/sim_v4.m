%%Ev
%%no test=0
%%ab=1
%%small=2
%%large=3
%%HCC
%%ab=1
%%pre=2
%%bleeding=1
%%death=4
%%transp=1
%%resect=2
%Need function for random variable HCC(seed,phcc) and EV(seed,pev) 
%and the output for these two random variables should contain seed and the
%states value described as above
%surgery,1-transplant, 2-resect
% GI bleeding after NSBB (20%) and banding ligation (14%). We could take the average for both (17%). 
% 
% Reference: Gluud LL, Krag A. Banding ligation versus beta-blockers for primary prevention in oesophageal varices in adults. Cochrane Database Syst Rev 2012:CD004544.



function [seed,cost,cost_egd_t,cost_hcc,cost_bb,cost_ebl,cost_trans,...
    cost_rsct,death_ev,death_hcc,death_cp,bleeding,surgery]=sim_v4(seed,...
    condition,Tmax,c_bb, c_hcc, p_bb,c_egd, c_ebl,c_trans, c_rsct,c_MRE,...
    speHcc, senHcc,hcc_inc,r_bld_suto,r_bld_luto,r_bld_lbbo,r_bld_leblo,...
mor_bld,mor_hcco,mor_cp,p_ev,p_evs,p_evto,...
p_hccpreo)
% CLDQ_non=1;
% CLDQ=0;
cycle=6;
clock=0;
bleeding=0;
death_ev=0;
death_hcc=0;
death_cp=0;
cost_egd=c_egd;
cost_hcc=c_hcc;
cost_egd2=0;
cost_bb=0;
cost_ebl=0;
cost_trans=0;
cost_rsct=0;
cost_conf=0;
% QALY=0;


inflation=1.03;
surgery=zeros(2);
state=zeros(Tmax,4);%ev,hcc,bleeding, treatment,transplant/resect
hidden_hcc=zeros(Tmax,1);
[seed, ini_state]=test_v2(seed,'EV',0,clock,hcc_inc,condition,p_ev,p_evs,p_evto,...
p_hccpreo);%EV
state(1,1)=ini_state;
[seed, ini_state]=test_v2(seed,'HCC',0,clock,hcc_inc,condition,p_ev,p_evs,p_evto,...
p_hccpreo);%HCC
hidden_hcc(1)=ini_state;

%determine obs of hcc
p_obsab_ab=speHcc/(speHcc+1-senHcc);
p_obspre_pre=senHcc/(1-speHcc+senHcc);
[seed,u]=u16807d(seed); 
if hidden_hcc(1)==1&&u<=p_obsab_ab
    state(1,2)=1;
elseif hidden_hcc(1)==1&&u>p_obsab_ab
    state(1,2)=1;
    cost_conf=cost_conf+c_MRE;
elseif hidden_hcc(1)==2&&u>p_obspre_pre
    state(1,2)=1;
else
    state(1,2)=2;
    cost_conf=cost_conf+c_MRE;
end
ti=2;
prior_ev=[state(1,1) 0];%ev state, time, bleeding, (0,1bb,2ebl)
prior_hcc=[state(1,2) 0];%1-ab,2-pre


%determine wether or not bb
[seed,u]=u16807d(seed); 
if u<p_bb
treat=1;
else
treat=2;
end
%Death transition probability (cycle 6 months)
r_bld_sut=1-exp(-ratecon(r_bld_suto(1),r_bld_suto(2))*cycle);
r_bld_lut=1-exp(-ratecon(r_bld_luto(1),r_bld_luto(2))*cycle);
r_bld_lbb=1-exp(-ratecon(r_bld_lbbo(1),r_bld_lbbo(2))*cycle);
r_bld_lebl=1-exp(-ratecon(r_bld_leblo(1),r_bld_leblo(2))*cycle);
mor_bld_p=1-exp(-ratecon(mor_bld(1),mor_bld(2))*cycle);
mor_hcc=1-exp(-ratecon(mor_hcco(1),mor_hcco(2))*cycle);
mor_cp_p=1-exp(-ratecon(mor_cp(1),mor_cp(2))*cycle);
%6*Tmax
while clock<Tmax

    clock=clock+cycle;
%     if condition==1
% %         QALY = QALY+cycle/12*qly_c;
%         if state(ti-1,1)>1 && state(ti-1,2)<2           
%             CLDQ = CLDQ + cycle/12*CLDQ_d;
%         elseif state(ti-1,2)>1
%             CLDQ = CLDQ + cycle/12*CLDQ_hcc;
%         elseif state(ti-1,1)==1 && state(ti-1,2)==1
%             CLDQ = CLDQ + cycle/12*CLDQ_c;
%         end
%     else
% %         QALY = QALY+cycle/12*1;
%         CLDQ = CLDQ + cycle/12*CLDQ_non;
%     end
    %t=clock;
    
    %mor_trans=1-exp(-ratecon(mor_transo(1),mor_transo(2))*cycle);
   if surgery(1)==0  %people without transplant have to be checked
       if (prior_ev(1)==1&&clock-prior_ev(2)==36) %EV test check
          %[seed,u]=u16389(seed); 
          if condition==1       
             [seed, d_state]=test_v2(seed,'EV',prior_ev(1),hcc_inc,36,condition,p_ev,p_evs,p_evto,...
    p_hccpreo);%EV test  
             state(ti,1)=d_state;
             prior_ev(1)=d_state;
             prior_ev(2)=clock;
          else
             state(ti,1)=1;
             prior_ev(1)=1;
             prior_ev(2)=clock;
          end
             cost_egd=cost_egd+c_egd/inflation^ceil(clock/12);
       elseif (prior_ev(1)==2&&clock-prior_ev(2)==24)
         [seed, d_state]=test_v2(seed,'EV',prior_ev(1),hcc_inc,24,condition,p_ev,p_evs,p_evto,...
p_hccpreo);%EV test  
         state(ti,1)=d_state;
         prior_ev(1)=d_state;
         prior_ev(2)=clock;
         cost_egd=cost_egd+c_egd/inflation^ceil(clock/12);

       elseif prior_ev(1)==3
           state(ti,1)=3;
           prior_ev(1)=3;
           prior_ev(2)=clock;
           cost_egd2=cost_egd2+c_egd/inflation^ceil(clock/12);
             if treat==1
                 state(ti-1,4)=1;
                 cost_bb=cost_bb+c_bb*cycle/inflation^ceil(clock/12);
             else
                 state(ti-1,4)=2;
                 cost_ebl=2.5*c_ebl/inflation^ceil(clock/12);
             end
       else
           state(ti,1)=state(ti-1,1);
       end
    %bleeding
    if condition==1
       %small bleeding
           if prior_ev(1)==2
             [seed,u]=u16807d(seed);
             if u<=r_bld_sut
                state(ti-1,3)=1;
                bleeding=1;
             elseif ti>2&&state(ti-2,3)==1
                state(ti-1,3)=1;
                bleeding=1;
             end  
           elseif prior_ev(1)==3
           %large bleeding
             [seed,u]=u16807d(seed);
                 if (state(ti-1,4)==0 && u<=r_bld_lut)||(state(ti-1,4)==1 && u<=r_bld_lbb)...
                         ||(state(ti-1,4)==2 && u<=r_bld_lebl)
                    state(ti-1,3)=1;
                    bleeding=1;
                 elseif ti>2&&state(ti-2,3)==1
                     state(ti-1,3)=1;
                     bleeding=1;
                 else
                 state(ti-1,3)=0;
                 bleeding=0;
                 end
           end
    end
%      if surgery(2)==0
       if prior_hcc(1)==1%HCC test check
           if condition ==1
              [seed, hidden_hcc(ti)]=test_v2(seed,'HCC',hidden_hcc(ti-1),hcc_inc,6,condition,p_ev,p_evs,p_evto,...
    p_hccpreo);%HCC test  
             [seed,u]=u16807d(seed); 
                if hidden_hcc(ti)==1&&u<=p_obsab_ab
                    state(ti,2)=1;
                elseif hidden_hcc(ti)==1&&u>p_obsab_ab
                    state(ti,2)=1;
                    cost_conf=cost_conf+c_MRE/inflation^ceil(clock/12);
                elseif hidden_hcc(ti)==2&&u>p_obspre_pre
                    state(ti,2)=1;
                else
                    state(ti,2)=2;
                    cost_conf=cost_conf+c_MRE/inflation^ceil(clock/12);
                end
              prior_hcc(1)=state(ti,2);
              prior_hcc(2)=clock;
           else
              state(ti,2)=1;
              prior_hcc(1)=state(ti,2);
              prior_hcc(2)=clock;
           end
          cost_hcc=cost_hcc+c_hcc/inflation^ceil(clock/12); 
       elseif (prior_ev(1)==1&&state(ti-1,2)==2)
%            if surgery(2)==0
                surgery(2)=1;
                cost_rsct=cost_rsct+c_rsct/inflation^ceil(clock/12);
                state(ti,2)=1;
                prior_hcc(1)=1;
                hidden_hcc(ti)=1;
%            end     
       elseif (prior_ev(1)==2&&state(ti-1,2)==2)
           surgery(1)=1;
           cost_trans=c_trans/inflation^ceil(clock/12);
           bleeding=0;
            break
       elseif (prior_ev(1)==3&&state(ti-1,2)==2)    
           surgery(1)=1;
           cost_trans=c_trans/inflation^ceil(clock/12);
           bleeding=0;
            break
       end
  
   else
       state(ti,1)=1;
       state(ti,2)=1;
   end  
  
    %Death
    if condition==1
        if hidden_hcc(ti-1)==2
            if surgery(1)==0
              if  state(ti-1,3)==1 %bleeding caused& hcc
        %           if condition==1% hidden_hcc(ti-1)==1 bleeding 
                     [seed,u]=u16807d(seed);
                      if u<=0.5
                          [seed,u]=u16807d(seed);
                          if u<= mor_bld_p
                          death_ev=1;
                          break 
                          end                     
                      else%hcc caused
                             [seed,u]=u16807d(seed);
                             if u<=mor_hcc
                                 death_hcc=1;
                                 break
                             end
                      end
              else
                     [seed,u]=u16807d(seed);
                     if u<=mor_hcc
                         death_hcc=1;
                         break
                     end
               end
            end
        elseif state(ti-1,3)==1 && hidden_hcc(ti-1)==1%bleeding 
            if state(ti-1,1)==2||state(ti-1,1)==3%bleeding  with EV small or large
              [seed,u]=u16807d(seed); 
              if u<=mor_bld_p
                  death_ev=1;
              break 
              end
            end
        elseif hidden_hcc(ti-1)==1 && state(ti-1,1)==1
            [seed,u]=u16807d(seed); 
            if u<=mor_cp_p
                death_cp=1;
                break
            end
        end     
    end
     
    ti=ti+1;    
end
% if surgery(1)==1
%     if Tmax-clock>12
%         CLDQ=CLDQ+CLDQ_trans+CLDQ_non*(Tmax-clock-12)/12;
%     else
%         CLDQ=CLDQ+CLDQ_trans*(Tmax-clock)/12;
%     end
% elseif surgery(2)
%     CLDQ=CLDQ+CLDQ_non*(Tmax-clock)/12;
% end
cost_egd_t=cost_egd+cost_egd2;
cost_hcc=cost_hcc+cost_conf;
cost=cost_egd+cost_egd2+cost_hcc+cost_bb+cost_ebl+cost_trans+cost_rsct+cost_conf;

end