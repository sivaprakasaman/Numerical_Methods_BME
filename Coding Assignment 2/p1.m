%Andrew Sivaprakasam
%BME 695 | Numerical Methods
%Coding Assignment 2 Problem 1

%Verification of spline functions

x1 = [0,1:3];
y1 = [0,3,4,7];
x2 = [2,3,5,6,9];
y2 = [4,6,7,5,-2];

figure;
title('Example A: n = 3');
hold on;
x_out = 1:0.01:3;
scatter(x1,y1);

[x_out,y_out] = spline_2(x1,y1,x_out,'natural');
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x1,y1,x_out,'complete');
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x1,y1,x_out,'parabolic');
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x1,y1,x_out,'nak');
plot(x_out,y_out,'linewidth',1.5);

legend('points','natural','complete','parabolic','nak','location','northwest');

hold off;
grid on
xlim([1,3]);

figure;
title('Example B: n = 5');
hold on;
x_out = 2:0.01:9;
scatter(x2,y2);

[x_out,y_out] = spline_2(x2,y2,x_out,'natural');
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x2,y2,x_out,'parabolic');
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x2,y2,x_out,'nak');
plot(x_out,y_out,'linewidth',1.5);

legend('points','natural','parabolic','nak');

hold off;
grid on

%% Runge Function

clear all
%generate sampled function
x = -1:0.2:1;
y = 1./(1+25*x.^2);
x_cont = -1:0.01:1;
y_cont = 1./(1+25*x_cont.^2);
x_out = x_cont;

figure;
hold on;
scatter(x,y);
plot(x_cont,y_cont,'k','linewidth',2.5);

[x_out,y_out] = spline_2(x,y,x_out,'natural')
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x,y,x_out,'complete')
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x,y,x_out,'parabolic')
plot(x_out,y_out,'linewidth',1.5);

[x_out,y_out] = spline_2(x,y,x_out,'nak')
plot(x_out,y_out,'linewidth',1.5);

P = polyfit(x,y,14);
y_out = polyval(P,x_out);
plot(x_out,y_out,'linewidth',1.5);

hold off

legend('sampled','continuous','natural','complete','parabolic','nak','polynomial_{14}');

title('Runge Function Comparisons of End Conditions');
xlim([-.7,.7])
grid on