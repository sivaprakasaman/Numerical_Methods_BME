function err = simulation_error(obs,params,alph)

l1 = num2str(params(1));
l2 = num2str(params(2));
l3 = num2str(params(3));

outputTemp = python('main.py',l1,l2,l3); %e.g. l1:'1.5', l2:'-2.2', l3:'1' 
output = sscanf(outputTemp,'%f'); %

H = output(1:100);
ICU = output(101:200);

mse_H = mean((H-obs(:,1)).^2);
mse_ICU = mean((ICU-obs(:,2)).^2);

%alpha weights the H error more. 
err = alph*mse_H + (1-alph)*mse_ICU; 

disp('iteration complete')

end

