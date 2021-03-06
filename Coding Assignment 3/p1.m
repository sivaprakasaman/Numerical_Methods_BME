%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 3

%% Implement Euler's to solve Lotka Volterra Model:

%essentially derive an exact solution to the model based on previous
%iteration and hold that for the time dt. 

%helpful resource! https://services.math.duke.edu/education/ccp/materials/diffcalc/predprey/pred3.html

%obviously this estimate will be poor if dt is large!!

b_1 = 2; % prey birthrate
p_1 = 0.04; % prob that preys pop size will be reduced by predator
p_2 = 0.02; % prob that predator pop size will grow by prey
m_2 = 1.06; % mortality rate
N1_0 = 100; % prey pop size
N2_0 = 15;  % pred pop size

% Time parameters (years)
dt = 0.01;
tmin = 0; 
tmax = 20;

t = (tmin + dt):dt:tmax;

N1 = zeros(1,length(t));
N2 = zeros(1,length(t));

N1(1) = N1_0;
N2(1) = N2_0;

for i = 2:length(t)
    
    N1(i) = N1(i-1) + (b_1*N1(i-1)-p_1*N1(i-1)*(N2(i-1)))*dt;
    N2(i) = N2(i-1) + (p_2*N1(i-1)*(N2(i-1))-m_2*N2(i-1))*dt;

end

hold on 
plot(t,N1,'linewidth',1.5);
plot(t,N2,'linewidth',1.5);
title('Uncalibrated Lotka Volterra Model');
xlabel('Time (yrs)');
ylabel('Population');
legend('N1','N2');

%% Model Calibration:

%Simple calibration...since N1/N2 are fixed, vary p_1 and p_2 such that min
%square error is acheived.

params0 = [b_1,p_1,p_2,m_2];
obs_t = [17.2];
obs_N12 = [29.7;51.8];
inits = [N1_0,N2_0];

params = fminsearch(@(params)LVerr(obs_t, obs_N12, params, inits),params0);

%re-run Euler's with the matched parameters:
for i = 2:length(t)
    
    N1(i) = N1(i-1) + (params(1)*N1(i-1)-params(2)*N1(i-1)*(N2(i-1)))*dt;
    N2(i) = N2(i-1) + (params(3)*N1(i-1)*(N2(i-1))-params(4)*N2(i-1))*dt;

end

figure;
hold on 
h1 = plot(t,N1,'linewidth',1.5);
datatip(h1, obs_t,obs_N12(1));
h2 = plot(t,N2,'linewidth',1.5);
datatip(h2,obs_t,obs_N12(2));
legend('N1','N2');
title('Calibrated Lotka Volterra Model');
xlabel('Time (yrs)');
ylabel('Population');
hold off