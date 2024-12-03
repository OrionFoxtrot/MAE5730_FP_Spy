clear; clc;

syms x1 y1 theta1 x2 y2 theta2 x3 y3 theta3 x4 y4 real
syms x1_dot y1_dot theta1_dot x2_dot y2_dot theta2_dot x3_dot y3_dot theta3_dot x4_dot y4_dot real
syms theta1_ddot x2_ddot y2_ddot theta2_ddot x3_ddot y3_ddot theta3_ddot x4_ddot y4_ddot real
syms rx1 ry1 rx2 ry2 rx3 ry3 d1 d2 d3 g m1 m2 m3 m4 ig1 ig2 ig3 real
syms F x1_ddot y1_ddot real

%Constraints
c1 = x1_ddot - d1*sin(theta1)*theta1_dot^2 + d1*cos(theta1)*theta1_ddot == x2_ddot; %x1 
c2 = d1*cos(theta1)*theta1_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot == y2_ddot;
c3 = x1_ddot - d1*sin(theta1)*theta1_dot^2 - d2*sin(theta2)*theta2_dot^2 + d1*cos(theta1)*theta1_ddot + d2*cos(theta2)*theta2_ddot == x3_ddot;
c4 = d1*cos(theta1)*theta1_dot^2 + d2*cos(theta2)*theta2_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot + d2*sin(theta2)*theta2_ddot == y3_ddot;
c5 = x1_ddot - d1*sin(theta1)*theta1_dot^2 - d2*sin(theta2)*theta2_dot^2 - d3*sin(theta3)*theta3_dot^2 + d1*cos(theta1)*theta1_ddot + d2*cos(theta2)*theta2_ddot + d3*cos(theta3)*theta3_ddot == x4_ddot;
c6 = d1*cos(theta1)*theta1_dot^2 + d2*cos(theta2)*theta2_dot^2 + d3*cos(theta3)*theta3_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot + d2*sin(theta2)*theta2_ddot + d3*sin(theta3)*theta3_ddot == y4_ddot;
c7 = x1_ddot == y1_ddot;

%NE:
n0 = rx1+F == m1*x1_ddot;
n00 = -ry1 - m1*g == m1*y1_ddot;

n1 = -rx1+rx2 == m2*x2_ddot;
n2 = ry1-ry2 -m2*g == m2*y2_ddot;
n3 = rx1*d1*cos(theta1)-ry1*d1*sin(theta1)+rx2*d1*cos(theta1)-ry2*d1*sin(theta1)== ig1 * theta1_ddot;

n4 = -rx2+rx3 == m3*x3_ddot;
n5 = ry2-ry3-m3*g==m3*y3_ddot;
n6 = rx2*d2*cos(theta2)-ry2*d2*sin(theta2)+rx3*d2*cos(theta2)-ry2*d2*sin(theta2)== ig2*theta2_ddot;

n7 = -rx3 == m3*x4_ddot;
n8 = ry3-m4*g == m4*y4_ddot;
n9 = rx3*d3*cos(theta3)-ry3*d3*sin(theta3) ==  ig3*theta3_ddot;

vars = [x1_ddot y1_ddot...
    x2_ddot y2_ddot theta1_ddot ...
    x3_ddot y3_ddot theta2_ddot ...
    x4_ddot y4_ddot theta3_ddot ...
    rx1 ry1 rx2 ...
    ry2 rx3 ry3 ];
eqns = [n0 n00 n1 n2 n3 n4 n5 n6 n7 n8 n9 c1 c2 c3 c4 c5 c6 c7];
[A,b] = equationsToMatrix(eqns,vars)

rank(A)
rank(b)
matlabFunction(A,'File','DAE_A_funs')
matlabFunction(b,'File','DAE_B_funs')


%%

clear; clc;

m1 = 1;
m2 = 1;
m3 = 1;
m4 = 1;
d1 = 5;
d2 = 5;
d3 = 5;
G = 9.81;

ig1 = m2*(d1*2)^2/12;
ig2 = m3*(d2*2)^2/12;
ig3 = m4*(d3*2)^2/12;

ff = pi/4;

X0_0 = [0,0]; % x1 y1
X0_1 = [0,0,0]; %x2 y2 theta1
X0_2 = [0,0,0]; %x3 y3 theta2
X0_3 = [0,0,0]; %x4 y4 theta3

X0_4 = [0,0]; % x1_dot y1_dot
X0_5 = [0,0,0]; %x2 y2 theta1 _DOTS
X0_6 = [0,0,0]; %x3 y3 theta2 _DOTS
X0_7 = [0,0,0]; %x4 y4 theta3 _DOTS

F = 0;

X0 = [X0_0, X0_1, X0_2, X0_3, X0_4, X0_5, X0_6, X0_7];

% Create time span for use with ode45

x_ddot_log = [];
y_ddot_log = [];

t = [0 15];
% Create function f(t,X) for use with ode45
fdynamic    = @(t,X) spy(t,X,m1,m2,m3,m4, d1/2,d2/2,d3/3,G,ig1,ig2,ig3, F);
% Solve the EoM with ode45
[t,X] = ode45(fdynamic,t,X0);


function Xdot = spy(t,X, m1,m2,m3,m4, d1,d2,d3,G, ig1,ig2,ig3, F)

x1 = X(1);
y1 = X(2);

x2 = X(3);
y2 = X(4);
theta1 = X(5);

x3 = X(6);
y3 = X(7);
theta2 = X(8);

x4 = X(9);
y4 = X(10);
theta3 = X(11);

x1_dot = X(12);
y1_dot = X(13);

x2_dot = X(14);
y2_dot = X(15);
theta1_dot = X(16);

x3_dot = X(17);
y3_dot = X(18);
theta2_dot = X(19);

x4_dot = X(20);
y4_dot = X(21);
theta3_dot = X(22);

%A = DAE_A_funs(D1,D2,D3,IG1,IG2,IG3,M1,M2,M3,M4,THETA1,THETA2,THETA3)
A = DAE_A_funs(d1,d2,d3,ig1,ig2,ig3,m1,m2,m3,m4,theta1,theta2,theta3);
%B = DAE_B_funs(F,D1,D2,D3,G,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT)
b = DAE_B_funs(F,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);

u = A\b;

x1_ddot = u(1)
y1_ddot = u(2)

x2_ddot = u(3);
y2_ddot = u(4);
theta1_ddot = u(5);

x3_ddot = u(6);
y3_ddot = u(7);
theta2_ddot = u(8);

x4_ddot = u(9);
y4_ddot = u(10);
theta3_ddot = u(11);

Xdot = [x1_dot,y1_dot, x2_dot y2_dot theta1_dot x3_dot y3_dot theta2_dot x4_dot y4_dot theta3_dot, ...
    x1_ddot,y1_ddot, x2_ddot y2_ddot, theta1_ddot, x3_ddot, y3_ddot, theta2_ddot, x4_ddot, y4_ddot, theta3_ddot]';
end