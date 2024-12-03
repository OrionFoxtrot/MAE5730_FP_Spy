
clear
clc
% L
syms xg yg theta1 theta2 theta3 d1 d2 d3 m1 m2 m3 m4 g real
syms xg_dot yg_dot theta1_dot theta2_dot theta3_dot real
syms xg_ddot yg_ddot theta1_ddot theta2_ddot theta3_ddot real
syms lam real
% xg_dot = diff(xg);
% yg_dot = diff(yg);
% theta1_dot = diff(theta1);
% theta2_dot = diff(theta2);
% theta3_dot = diff(theta3);
 

Ig1 = 1/12*m2*(2*d1)^2;
Ig2 = 1/12*m3*(2*d2)^2;
Ig3 = 1/12*m4*(2*d3)^2;

vg = sqrt(xg_dot^2+yg_dot^2);

vg2 = sqrt((xg_dot+d1*theta1_dot*cos(theta1))^2 + (yg_dot-d1*theta1_dot*sin(theta1))^2);
vg3 = sqrt((xg_dot+d1*theta1_dot*cos(theta1)+d2*theta2_dot*cos(theta2))^2 + (yg_dot-d1*theta1_dot*sin(theta1)-d2*theta2_dot*sin(theta2))^2);
vg4 = sqrt((xg_dot+d1*theta1_dot*cos(theta1)+d2*theta2_dot*cos(theta2)+d3*theta3_dot*cos(theta3))^2 + (yg_dot-d1*theta1_dot*sin(theta1)-d2*theta2_dot*sin(theta2)-d3*theta3_dot*sin(theta3))^2);

% KE = 1/2* ( m1*vg^2 ...
%     +m2*d1^2*theta1_dot^2+Ig1*d1^2*theta1_dot^2 ...
%     +m3*d2^2*theta2_dot^2+Ig2*d2^2*theta2_dot^2 ...
%     +m4*d3^2*theta3_dot^2+Ig3*d3^2*theta3_dot^2);
KE = 1/2* ( m1*vg^2 ...
    +m2*vg2^2+Ig1*theta1_dot^2 ...
    +m3*vg3^2+Ig2*theta2_dot^2 ...
    +m4*vg4^2+Ig3*theta3_dot^2 );
PE = g* (m1*yg ...
    + m2*(yg-d1*cos(theta1))...
    + m3*(yg-d1*cos(theta1)-d2*cos(theta2))...
    + m4*(yg - d1*cos(theta1)-d2*cos(theta2)-d3*cos(theta3)));

L = KE-PE;


J = @(f,x) jacobian(f,x);
q = [xg yg theta1 theta2 theta3]';
q_dot = [xg_dot yg_dot theta1_dot theta2_dot theta3_dot]';
q_ddot = [xg_ddot yg_ddot theta1_ddot theta2_ddot theta3_ddot]';

Q = [100-lam, lam, -0.2, -0.2, -0.2]'; %100 is forcing lam is constraint

%Q = [0,0,0,0,0]';

%EoM = J(J(L,q_dot),q)*q_dot + J(J(L,q_dot),q_dot)*q_ddot - J(L,q) == Q;
EoM = J(J(L,q_dot),q_dot)*q_ddot +  J(J(L,q_dot),q)*q_dot - J(L,q)' == Q;

eom1 = EoM(1);
eom2 = EoM(2);
eom3 = EoM(3);
eom4 = EoM(4);
eom5 = EoM(5);



sol = solve([eom1;eom2;eom3;eom4;eom5],[xg_ddot,yg_ddot,theta1_ddot,theta2_ddot,theta3_ddot])

%constraint 1
%y_ddot-x_ddot = 0
lamfun = sol.yg_ddot-sol.xg_ddot == 0;
lam_sol = solve(lamfun,lam);

matlabFunction([sol.xg_ddot sol.yg_ddot sol.theta1_ddot, sol.theta2_ddot,sol.theta3_ddot],'File','sol_funs')
matlabFunction(lam_sol,'File','lam_func')

%%
clear; clc;

m1 = 1;
m2 = 5;
m3 = 1;
m4 = 1;
d1 = 5;
d2 = 5;
d3 = 5;
G = 9.81;

ff = pi/4;
X0_0 = [0,0,0,0,0]; %xg yg theta1 theta2 theta3
X0_1 = [0,0,0,0,0]; %xg_dot yg_dot, theta1_dot, theta2_dot, theta3_dot
X0 = [X0_0, X0_1];

% Create time span for use with ode45


t = [0 15];
% Create function f(t,X) for use with ode45
fdynamic    = @(t,X) spy(t,X,m1,m2,m3,m4, d1,d2,d3,G);
% Solve the EoM with ode45
[t,X] = ode45(fdynamic,t,X0);


function Xdot = spy(t,X, m1,m2,m3,m4, d1,d2,d3,G)

xg = X(1);
yg = X(2);
theta1 = X(3);
theta2 = X(4);
theta3 = X(5);
xg_dot = X(6);
yg_dot = X(7);
theta1_dot = X(8);
theta2_dot = X(9);
theta3_dot = X(10);

%LAM_SOL = LAM_SOL = LAM_FUNC(D1,D2,D3,G,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT)
lam =      lam_func(d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);
%OUT1 = SOL_FUNS(D1,D2,D3,G,LAM,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT)
ddots = sol_funs(d1,d2,d3,G,lam,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);

Xdot = [xg_dot, yg_dot, theta1_dot, theta2_dot, theta3_dot, ddots]';
end