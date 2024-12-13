%% Equations of Motion Generation

clear; clc;
function gen = generate_matrix(stage)
    syms x1 y1 theta1 x2 y2 theta2 x3 y3 theta3 x4 y4 real
    syms x1_dot y1_dot theta1_dot x2_dot y2_dot theta2_dot x3_dot y3_dot theta3_dot x4_dot y4_dot real
    syms theta1_ddot x2_ddot y2_ddot theta2_ddot x3_ddot y3_ddot theta3_ddot x4_ddot y4_ddot real
    syms rx1 ry1 rx2 ry2 rx3 ry3 d1 d2 d3 g m1 m2 m3 m4 ig1 ig2 ig3 real
    syms F F_N_x F_N_y x1_ddot y1_ddot c real
    
    %Constraints
    % c1 = x1_ddot - d1*sin(theta1)*theta1_dot^2 + d1*cos(theta1)*theta1_ddot == x2_ddot; %x1
    % c2 = d1*cos(theta1)*theta1_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot == y2_ddot;
    % c3 = x1_ddot - d1*sin(theta1)*theta1_dot^2 - d2*sin(theta2)*theta2_dot^2 + d1*cos(theta1)*theta1_ddot + d2*cos(theta2)*theta2_ddot == x3_ddot;
    % c4 = d1*cos(theta1)*theta1_dot^2 + d2*cos(theta2)*theta2_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot + d2*sin(theta2)*theta2_ddot == y3_ddot;
    % c5 = x1_ddot - d1*sin(theta1)*theta1_dot^2 - d2*sin(theta2)*theta2_dot^2 - d3*sin(theta3)*theta3_dot^2 + d1*cos(theta1)*theta1_ddot + d2*cos(theta2)*theta2_ddot + d3*cos(theta3)*theta3_ddot == x4_ddot;
    % c6 = d1*cos(theta1)*theta1_dot^2 + d2*cos(theta2)*theta2_dot^2 + d3*cos(theta3)*theta3_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot + d2*sin(theta2)*theta2_ddot + d3*sin(theta3)*theta3_ddot == y4_ddot;
    
    c1 = x1_ddot - d1*sin(theta1)*theta1_dot^2 + d1*cos(theta1)*theta1_ddot == x2_ddot;
    c2 = d1*cos(theta1)*theta1_dot^2 + y1_ddot + d1*sin(theta1)*theta1_ddot == y2_ddot;
    c3 = x1_ddot - 2*d1*sin(theta1)*theta1_dot^2 - d2*sin(theta2)*theta2_dot^2 + 2*d1*cos(theta1)*theta1_ddot + d2*cos(theta2)*theta2_ddot == x3_ddot;
    c4 = 2*d1*cos(theta1)*theta1_dot^2 + d2*cos(theta2)*theta2_dot^2 + y1_ddot + 2*d1*sin(theta1)*theta1_ddot + d2*sin(theta2)*theta2_ddot == y3_ddot;
    c5 = x1_ddot - 2*d1*sin(theta1)*theta1_dot^2 - 2*d2*sin(theta2)*theta2_dot^2 - d3*sin(theta3)*theta3_dot^2 + 2*d1*cos(theta1)*theta1_ddot + 2*d2*cos(theta2)*theta2_ddot + d3*cos(theta3)*theta3_ddot == x4_ddot;
    c6 = 2*d1*cos(theta1)*theta1_dot^2 + 2*d2*cos(theta2)*theta2_dot^2 + d3*cos(theta3)*theta3_dot^2 + y1_ddot + 2*d1*sin(theta1)*theta1_ddot + 2*d2*sin(theta2)*theta2_ddot + d3*sin(theta3)*theta3_ddot == y4_ddot;
    
    
    % Func y = x
    if stage == 0
        c7 = x1_ddot == y1_ddot;
        c8 = F_N_x + F_N_y == 0;
    end
    
    % Func y = x^2
    if stage == 90
        c7 = 2*x1_dot^2 + 2*x1*x1_ddot == y1_ddot;
        c8 = F_N_x + 2*x1*F_N_y == 0;
    end
    
    % Func y = sin(x)
    if stage == 91
        c7 = y1_ddot == 2/3*cos(x1/3)*x1_ddot - 2/9*sin(x1/3)*(x1_dot)^2;
        c8 = F_N_x + 2/3*cos(x1/3)*F_N_y == 0;
    end
    
    if stage == 91
        c7 = -2/9*cos(x1/3)*x1_dot^2 - 2/3 * sin(x1/3)*x1_ddot == y1_ddot;
        c8 = F_N_x - 2/3*sin(x1/3)*F_N_y == 0;
    end
    
    % y = (x/10-3)^2 -1
    if stage == 0
        c7 = y1_ddot == x1_dot^2/50 + (x1/10 - 3)*x1_ddot/5;
        c8 = F_N_x+ (x1/50 - 3/5)*F_N_y == 0;
    end
    
    % y = -1/10x+16 [60,100]
    if stage == 1
        c7 = y1_ddot == -1/10*x1_ddot;
        c8 = F_N_x+ (-1/10)*F_N_y == 0;
    end
    if stage == 2
        c7 = y1_ddot == -(x1_dot)^2/50 - (x1/10-12)*x1_ddot/5;
        c8 = F_N_x + (12/5 - x1/50)*F_N_y == 0;
    end
    
    %NE:
    n0 = rx1+F + F_N_x == m1*x1_ddot;
    n00 = -ry1 - m1*g + F_N_y == m1*y1_ddot;
    
    n1 = -rx1+rx2 == m2*x2_ddot;
    n2 = ry1-ry2 -m2*g == m2*y2_ddot;
    n3 = rx1*d1*cos(theta1)-ry1*d1*sin(theta1)+rx2*d1*cos(theta1)-ry2*d1*sin(theta1) -c*theta1_dot == ig1 * theta1_ddot;
    
    n4 = -rx2+rx3 == m3*x3_ddot;
    n5 = ry2-ry3-m3*g==m3*y3_ddot;
    n6 = rx2*d2*cos(theta2)-ry2*d2*sin(theta2)+rx3*d2*cos(theta2)-ry3*d2*sin(theta2) -c*theta2_dot== ig2*theta2_ddot;
    
    n7 = -rx3 == m4*x4_ddot;
    n8 = ry3-m4*g == m4*y4_ddot;
    n9 = rx3*d3*cos(theta3)-ry3*d3*sin(theta3) -c*theta3_dot == ig3*theta3_ddot;
    
    vars = [x1_ddot y1_ddot x2_ddot y2_ddot theta1_ddot x3_ddot y3_ddot theta2_ddot x4_ddot y4_ddot theta3_ddot ...
        rx1 ry1 rx2 ry2 rx3 ry3 F_N_x F_N_y];
    eqns = [n0 n00 n1 n2 n3 n4 n5 n6 n7 n8 n9 c1 c2 c3 c4 c5 c6 c7 c8];
    [A,b] = equationsToMatrix(eqns,vars);
    
    rank(A)
    rank(b)
    v1 = sprintf('DAE_Matrix/DAE_A_funs_%s',num2str(stage));
    v2 = sprintf('DAE_Matrix/DAE_B_funs_%s',num2str(stage));

    matlabFunction(A,'File',v1);
    matlabFunction(b,'File',v2);

end

generate_matrix(2) % Generate Matrix for course section 3. (This uses standard programming indexing)
%% Complete Matrix Generation
generate_matrix(0);
generate_matrix(1);
generate_matrix(2);
disp('generation complete')
%% Simulation Block

clear; clc;
addpath("DAE_Matrix\")
m1 = 30;
m2 = 15;
m3 = 5;
m4 = 2;

d1 = 5;
d2 = 5;
d3 = 5;
G = 9.81;



ig1 = m2*(d1)^2/12;
ig2 = m3*(d2)^2/12;
ig3 = m4*(d3)^2/12;

ff = pi/4;

X0_0 = [0,10]; % x1 y1
X0_1 = [0,10-d1,0]; %x2 y2 theta1
X0_2 = [0,10-2*d1-d2,0]; %x3 y3 theta2
X0_3 = [0,10-2*d1-2*d2-d3,0]; %x4 y4 theta3

X0_4 = [0,0]; % x1_dot y1_dot
X0_5 = [0,0,0]; %x2 y2 theta1 _DOTS
X0_6 = [0,0,0]; %x3 y3 theta2 _DOTS
X0_7 = [0,0,0]; %x4 y4 theta3 _DOTS

c=1; %should be positive
F = 100;
X0 = [X0_0, X0_1, X0_2, X0_3, X0_4, X0_5, X0_6, X0_7];

t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy(t,X,m1,m2,m3,m4, d1/2,d2/2,d3/2,G,ig1,ig2,ig3, F,c);
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events',@stage0_stop);
[t0,X] = ode15s(fdynamic,t,X0, options);


lastX = X(height(X),:);
% lastX(12)=0;
% lastX(13)=0;
lastX(13) = stage2_dot(lastX(1),lastX(12));

F = 100;
t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy1(t,X,m1,m2,m3,m4, d1/2,d2/2,d3/2,G,ig1,ig2,ig3, F,c);
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events',@stage1_stop);
[t1,X1] = ode15s(fdynamic,t,lastX, options);

lastX1 = X1(height(X1),:);
% lastX1(12)=0;
% lastX1(13)=0;
lastX1(13) = stage3_dot(lastX1(1),lastX1(12));

F = 100;
t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy2(t,X,m1,m2,m3,m4, d1/2,d2/2,d3/2,G,ig1,ig2,ig3, F,c);
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events',@stage2_stop);
[t2,X2] = ode15s(fdynamic,t,lastX1, options);

X = [X;X1;X2];

function [value, isterminal, direction] = stage0_stop(t, X)
value      = (X(1) >= 60);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = stage1_stop(t, X)
value      = (X(1) >= 100);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = stage2_stop(t, X)
value      = (X(1) >= 140);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function Xdot = spy(t,X, m1,m2,m3,m4, d1,d2,d3,G, ig1,ig2,ig3, F,c)

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

%A = DAE_A_funs(D1,D2,D3,IG1,IG2,IG3,M1,M2,M3,M4,THETA1,THETA2,THETA3,X1)
A = DAE_A_funs_0(d1,d2,d3,ig1,ig2,ig3,m1,m2,m3,m4,theta1,theta2,theta3,x1);
%B = DAE_B_funs(F,D1,D2,D3,G,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT,X1,X1_DOT)
b = DAE_B_funs_0(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot, x1_dot);

u = A\b;

%vars = [x1_ddot y1_ddot x2_ddot y2_ddot theta1_ddot x3_ddot y3_ddot theta2_ddot x4_ddot y4_ddot theta3_ddot ...
%   rx1 ry1 rx2 ry2 rx3 ry3 F_N_x F_N_y];

x1_ddot = u(1);
y1_ddot = u(2);

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

function Xdot = spy1(t,X, m1,m2,m3,m4, d1,d2,d3,G, ig1,ig2,ig3, F,c)

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
A = DAE_A_funs_1(d1,d2,d3,ig1,ig2,ig3,m1,m2,m3,m4,theta1,theta2,theta3);
%B = DAE_B_funs(F,D1,D2,D3,G,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT)
b = DAE_B_funs_1(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);

u = A\b;

%vars = [x1_ddot y1_ddot x2_ddot y2_ddot theta1_ddot x3_ddot y3_ddot theta2_ddot x4_ddot y4_ddot theta3_ddot ...
%   rx1 ry1 rx2 ry2 rx3 ry3 F_N_x F_N_y];

x1_ddot = u(1);
y1_ddot = u(2);

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

function Xdot = spy2(t,X, m1,m2,m3,m4, d1,d2,d3,G, ig1,ig2,ig3, F,c)

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

%A = DAE_A_funs(D1,D2,D3,IG1,IG2,IG3,M1,M2,M3,M4,THETA1,THETA2,THETA3,X1)
A = DAE_A_funs_2(d1,d2,d3,ig1,ig2,ig3,m1,m2,m3,m4,theta1,theta2,theta3,x1);
%B = DAE_B_funs(F,D1,D2,D3,G,M1,M2,M3,M4,THETA1,THETA2,THETA3,THETA1_DOT,THETA2_DOT,THETA3_DOT,X1,X1_DOT)
b = DAE_B_funs_2(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot, x1_dot);

u = A\b;

%vars = [x1_ddot y1_ddot x2_ddot y2_ddot theta1_ddot x3_ddot y3_ddot theta2_ddot x4_ddot y4_ddot theta3_ddot ...
%   rx1 ry1 rx2 ry2 rx3 ry3 F_N_x F_N_y];

x1_ddot = u(1);
y1_ddot = u(2);

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

function y = stage1(x)
    y = (x/10-3)^2+1;
end
function y_dot = stage1_dot(x, x_dot)
    y_dot = -(x/10-12)*x_dot/5;
end

function y = stage2(x)
    y = -(1/10)*x+16;
end
function y_dot = stage2_dot(x, x_dot)
    y_dot = -x_dot/10;
end

function y = stage3(x)
    y = -(x/10-12)^2+10;
end
function y_dot = stage3_dot(x, x_dot)
    y_dot = -(x/10-12)*x_dot/5;
end


%%

% X0_0 = [1,2*cos(1/3)]; % x1 y1
% X0_1 = [1,X0_0(2)-d1/2,0]; %x2 y2 theta1
% X0_2 = [1,X0_0(2)-d1-d2/2,0]; %x3 y3 theta2
% X0_3 = [1,X0_0(2)-d1-d2-d3,0]; %x4 y4 theta3
