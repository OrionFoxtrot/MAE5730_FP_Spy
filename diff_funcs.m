clear; clc;
syms x1(t) y1(t) theta1(t) x2(t) y2(t) theta2(t) x3(t) y3(t) theta3(t) rx1 ry1 rx2 ry2 rx3 ry3 t
syms d1 d2 d3
syms x2_dot y2_dot x3_dot y3_dot x4_dot y4_dot
syms x2_ddot y2_ddot x3_ddot y3_ddot x4_ddot y4_ddot


x2 = d1*sin(theta1(t))+x1(t);
y2 = -d1*cos(theta1(t))+y1(t);

x3 = x2+d2*sin(theta2(t));
y3 = y2-d2*cos(theta2(t));

x4 = x3+d3*sin(theta3(t));
y4 = y3-d3*cos(theta3(t));
%%
clear; clc;
syms x1(t) y1(t) theta1(t) x2(t) y2(t) theta2(t) x3(t) y3(t) theta3(t) rx1 ry1 rx2 ry2 rx3 ry3 t
syms d1 d2 d3
syms x2_dot y2_dot x3_dot y3_dot x4_dot y4_dot
syms x2_ddot y2_ddot x3_ddot y3_ddot x4_ddot y4_ddot


x2 = x1(t) + d1*sin(theta1(t));
y2 = y1(t) - d1*cos(theta1(t));

x3 = x1(t) + 2*d1*sin(theta1(t)) + d2*sin(theta2(t));
y3 = y1(t) - 2*d1*cos(theta1(t)) - d2*cos(theta2(t));

x4 = x1(t) + 2*d1*sin(theta1(t)) + 2*d2*sin(theta2(t)) + d3*sin(theta3(t));
y4 = y1(t) - 2*d1*cos(theta1(t)) - 2*d2*cos(theta2(t)) - d3*cos(theta3(t));
%% DAE
clear; clc;
syms x(t) t

% y = x^2;
% pretty(simplify(diff(diff(y,t))))
% 
% y = 2*cos(x(t)/3);
% pretty(simplify(diff(diff(y,t))))
% pretty(simplify(diff(y,x)))
% 
% y = (x/10-3)^2-1
% pretty(simplify(diff(diff(y,t))))
% pretty(simplify(diff(y,x)))

y = -(x/10-12)^2+10
pretty(simplify(diff(diff(y,t))))
pretty(simplify(diff(y,x)))
%% EL
clear; clc;
syms x(t) t

y = -(x(t)/10-12)^2+10;
%y = cos(x(t));
pretty(simplify((diff(y,t))))
pretty(simplify(diff(diff(y,t))))

%%
clear; clc;

syms x1(t) y1(t) theta1(t) x2(t) y2(t) theta2(t) x3(t) y3(t) theta3(t) rx1 ry1 rx2 ry2 rx3 ry3 t
syms d1 d2 d3
syms x2_dot y2_dot x3_dot y3_dot x4_dot y4_dot
syms x2_ddot y2_ddot x3_ddot y3_ddot x4_ddot y4_ddot

x2 = x1(t) + d1*sin(theta1(t));
y2 = y1(t) - d1*cos(theta1(t));

x3 = x1(t) + 2*d1*sin(theta1(t)) + d2*sin(theta2(t));
y3 = y1(t) - 2*d1*cos(theta1(t)) - d2*cos(theta2(t));

x4 = x1(t) + 2*d1*sin(theta1(t)) + 2*d2*sin(theta2(t)) + d3*sin(theta3(t));
y4 = y1(t) - 2*d1*cos(theta1(t)) - 2*d2*cos(theta2(t)) - d3*cos(theta3(t));


%%
clear; clc;
syms c
vpasolve((-100/10+c)^2+1==6,c)