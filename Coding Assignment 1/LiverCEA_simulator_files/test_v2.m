%state ab=1
%state small=2
%state large=3
%state pre=2
function [seed,d_state]=test_v2(seed,test,state,hcc_inc,t,condition,p_ev,p_evs,p_evto,...
p_hccpreo)

% if ~exist('condition','var')
%     condition=true;
% end
           
[seed,u]=u16807d(seed); 

% p_ev=.653;
% p_evs=.48;

p_evt=1-exp(-p_evto*t);
p_hccpre=hcc_inc;
p_hccpret=1-exp(-p_hccpreo*t);


switch test
    case 'EV'  
        if condition==1
            if state==0
                if u>p_ev
                    d_state=1;
                else
                    [seed,u]=u16807d(seed); 
                    if u<=p_evs
                        d_state=2;
                    else
                        d_state=3;
                    end
                end
            elseif state==1
                 if u>p_evt
                    d_state=1;
                else
                    [seed,u]=u16807d(seed); 
                    if u<=p_evs
                        d_state=2;
                    else
                        d_state=3;
                    end
                end
            elseif state==2
                 if u<=p_evs
                    d_state=2;
                 else 
                    d_state=3;
                 end
            else
                 d_state=3;
            end
        else
            d_state=1;
        end
    case 'HCC'
        
        if condition==1
            if state==0
                if u>p_hccpre
                    d_state=1;
                else
                    d_state=2;
                end
            elseif state==1
             if u>p_hccpret
                d_state=1;
             else
                d_state=2;
             end
            else 
                d_state=2;
            end
        else 
            d_state=1;
        end
end   

end
