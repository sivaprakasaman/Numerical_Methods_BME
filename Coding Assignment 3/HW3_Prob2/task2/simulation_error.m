function err = simulation_error(obs,params,seed,alph)

l1 = num2str(params(1));
l2 = num2str(params(2));
l3 = num2str(params(3));

output = zeros(30,2);

parfor i = 1:30
    outputTemp = python('main.py',l1,l2,l3,num2str(seed(i))); %e.g. l1:'3.15', l2:'-0.2', l3:'8' 
    output(i,:) = transpose(sscanf(outputTemp,'%f')); %
end

H = output(:,1);
ICU = output(:,2);
% 
% mse_H = mean((H-obs(:,1)).^2);
% mse_ICU = mean((ICU-obs(:,2)).^2);
% 
% %alpha weights the H error more. 
% err = alph*mse_H + (1-alph)*mse_ICU; 


%Normality assumption tends to hold here, so maximizing p value of t test.
% 
% [~, p_H] = ttest2(H,obs(:,1));
% [~, p_ICU] = ttest2(ICU,obs(:,2));
% 
% 
% err = (alph)*(1-p_H) + (1-alph)*(1-p_ICU);s

mse_H = mean((H-obs(:,1)).^2);
mse_ICU = mean((ICU-obs(:,2)).^2);

%alpha weights the H error more. 
err = alph*mse_H + (1-alph)*mse_ICU; 

disp('iteration complete')

end

