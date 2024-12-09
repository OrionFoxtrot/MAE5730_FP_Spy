
clear
clc
% L
function gen = generate_equations(stage)
    syms xg yg theta1 theta2 theta3 d1 d2 d3 m1 m2 m3 m4 g real
    syms xg_dot yg_dot theta1_dot theta2_dot theta3_dot real
    syms xg_ddot yg_ddot theta1_ddot theta2_ddot theta3_ddot real
    syms lam c F real
    % xg_dot = diff(xg);
    % yg_dot = diff(yg);
    % theta1_dot = diff(theta1);
    % theta2_dot = diff(theta2);
    % theta3_dot = diff(theta3);
    
    
    Ig1 = 1/12*m2*(2*d1)^2;
    Ig2 = 1/12*m3*(2*d2)^2;
    Ig3 = 1/12*m4*(2*d3)^2;
    
    vg = sqrt(xg_dot^2+yg_dot^2);
    
    % vg2 = sqrt((xg_dot+d1*theta1_dot*cos(theta1))^2 + (yg_dot-d1*theta1_dot*sin(theta1))^2);
    % vg3 = sqrt((xg_dot+d1*2*theta1_dot*cos(theta1)+d2*theta2_dot*cos(theta2))^2 + (yg_dot-2*d1*theta1_dot*sin(theta1)-d2*theta2_dot*sin(theta2))^2);
    % vg4 = sqrt((xg_dot+d1*theta1_dot*cos(theta1)+d2*theta2_dot*cos(theta2)+d3*theta3_dot*cos(theta3))^2 + (yg_dot-d1*theta1_dot*sin(theta1)-d2*theta2_dot*sin(theta2)-d3*theta3_dot*sin(theta3))^2);
    
    vg2 = sqrt((xg_dot + d1*cos(theta1)*theta1_dot)^2 + (yg_dot + d1*sin(theta1)*theta1_dot)^2);
    vg3 = sqrt((xg_dot + 2*d1*cos(theta1)*theta1_dot + d2*cos(theta2)*theta2_dot)^2 + (yg_dot + 2*d1*sin(theta1)*theta1_dot + d2*sin(theta2)*theta2_dot)^2);
    vg4 = sqrt((xg_dot + 2*d1*cos(theta1)*theta1_dot + 2*d2*cos(theta2)*theta2_dot + d3*cos(theta3)*theta3_dot)^2 + (yg_dot + 2*d1*sin(theta1)*theta1_dot + 2*d2*sin(theta2)*theta2_dot + d3*sin(theta3)*theta3_dot)^2);
    
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
        + m3*(yg-2*d1*cos(theta1)-d2*cos(theta2))...
        + m4*(yg -2*d1*cos(theta1)-2*d2*cos(theta2)-d3*cos(theta3)));
    
    L = KE-PE;
    
    
    J = @(f,x) jacobian(f,x);
    q = [xg yg theta1 theta2 theta3]';
    q_dot = [xg_dot yg_dot theta1_dot theta2_dot theta3_dot]';
    q_ddot = [xg_ddot yg_ddot theta1_ddot theta2_ddot theta3_ddot]';
    
    % Q MUST BE REPLACED %100 is forcing lam is constraint
    %c = 1;
    
    if stage == 0
        Q = [F-1/5*(xg/10-3)*lam, lam, c*theta1_dot,c*theta1_dot,c*theta1_dot]'; % y = (x/10-3)^2+1
    end
    if stage == 1
        Q = [F-1/10*lam, lam, c*theta2_dot,c*theta2_dot, c*theta2_dot]'; % y = -1/10x+16
    end
    if stage == 2
        Q = [F+1/5*(xg/10-12)*lam, lam, c*theta3_dot, c*theta3_dot, c*theta3_dot]'; % y = -(x/10-12)^2+10
    end

    %250 aprox F
    
    %Q = [600-lam, lam, -0.2, -0.2, -0.2]'; % y = x
    %Q = [-(2*xg+3)*lam, lam,-0.2,-0.2,-0.2]'; % y = (x-3)^2+1
    %Q = [200-1/5*(xg/10-3)*lam, lam,-0.2,-0.2,-0.2]'; % y = (x/10-3)^2+1
    %Q = [1+sin(xg)*lam, lam, -0.2, -0.2,-0.2]'; % y = cos(x)
    %Q = [0,0,0,0,0]';
    
    EoM = J(J(L,q_dot),q_dot)*q_ddot +  J(J(L,q_dot),q)*q_dot - J(L,q)' == Q;
    
    eom1 = EoM(1);
    eom2 = EoM(2);
    eom3 = EoM(3);
    eom4 = EoM(4);
    eom5 = EoM(5);

    sol = solve([eom1;eom2;eom3;eom4;eom5],[xg_ddot,yg_ddot,theta1_ddot,theta2_ddot,theta3_ddot]);
    
    %constraint 1
    %y_ddot-x_ddot = 0 %LAMFUN MUST BE REPLACED WITH FUN
    
    if stage == 0
        lamfun = sol.yg_ddot == 1/50*xg_dot^2 + 1/5*(xg/10-3)*sol.xg_ddot; % y = (x/10-3)^2 + 1
    end
    if stage == 1
        lamfun = sol.yg_ddot+1/10*sol.xg_ddot == 0; % y = x
    end
    if stage == 2
        lamfun = sol.yg_ddot == -1/50*xg_dot^2 - 1/5*(xg/10-12)*sol.xg_ddot; %y = -(x/10-12)^2+10
    end
    %lamfun = sol.yg_ddot-sol.xg_ddot == 0; % y = x
    %lamfun = sol.yg_ddot - 2*xg_dot^2-2*(xg-3)*sol.xg_ddot == 0; % y = (x-3)^2 + 1
    %lamfun = sol.yg_ddot == 1/50*xg_dot^2 + 1/5*(xg/10-3)*sol.xg_ddot; % y = (x/10-3)^2 + 1
    %lamfun = sol.yg_ddot == -sin(xg)*sol.xg_ddot - cos(xg)*xg_dot^2 ; %y = cos(x)
    lam_sol = solve(lamfun,lam);
    
    v1 = sprintf('EL_Eqs/sol_funs_%s',num2str(stage));
    v2 = sprintf('EL_Eqs/lam_func_%s',num2str(stage));

    matlabFunction([sol.xg_ddot sol.yg_ddot sol.theta1_ddot, sol.theta2_ddot,sol.theta3_ddot],'File',v1)
    matlabFunction(lam_sol,'File',v2)
end

%single generation
generate_equations(0)

%% Course Generation
clear; clc;
generate_equations(0)
generate_equations(1)
generate_equations(2)
disp('generation complete')
%%
clear; clc;
addpath("EL_Eqs\");

m1 = 30;
m2 = 15;
m3 = 5;
m4 = 2;

d1 = 5;
d2 = 5;
d3 = 5;

G = 9.81;

c=-1; %Should be Neg

ff = pi/4;
X0_0 = [0, 10, 0,0,0];
X0_1 = [0,0,0,0,0]; %xg_dot yg_dot, theta1_dot, theta2_dot, theta3_dot
X0 = [X0_0, X0_1];


%Stage 0 => y = (x/10-3)^2+1
F = 100;
t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy(t,X,m1,m2,m3,m4,d1/2,d2/2,d3/2,G,c,F);
options    = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events', @stage0_stop);
[t,X] = ode45(fdynamic,t,X0,options);

lastX = X(height(X),:);
% lastX(6)=0;
% lastX(7)=0;
lastX(7)=stage2_dot(lastX(1), lastX(6));

%Stage 1 => y = -(1/10)x+16
F = 100;
t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy1(t,X,m1,m2,m3,m4,d1/2,d2/2,d3/2,G,c,F);
options    = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events', @stage1_stop);
[t1,X1] = ode45(fdynamic,t,lastX,options);

lastX1 = X1(height(X1),:);
% lastX1(6)=0;
% lastX1(7)=0;
lastX1(7)=stage3_dot(lastX1(1), lastX1(6));

%Stage 2 => -(x/10-12)^2+16
F = 100;
t = linspace(0,100,9e2);
fdynamic    = @(t,X) spy2(t,X,m1,m2,m3,m4,d1/2,d2/2,d3/2,G,c,F);
options    = odeset('RelTol',1e-9,'AbsTol',1e-9,'Refine',2,'Events', @stage2_stop);
[t2,X2] = ode45(fdynamic,t,lastX1,options);


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


function Xdot = spy(t,X, m1,m2,m3,m4, d1,d2,d3,G,c,F)

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
lam =      lam_func_0(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot,xg,xg_dot);
ddots = sol_funs_0(F,c,d1,d2,d3,G,lam,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot,xg);

Xdot = [xg_dot, yg_dot, theta1_dot, theta2_dot, theta3_dot, ddots]';
end

function Xdot = spy1(t,X, m1,m2,m3,m4, d1,d2,d3,G,c,F)

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

lam = lam_func_1(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);
ddots = sol_funs_1(F,c,d1,d2,d3,G,lam,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot);

Xdot = [xg_dot, yg_dot, theta1_dot, theta2_dot, theta3_dot, ddots]';
end

function Xdot = spy2(t,X, m1,m2,m3,m4, d1,d2,d3,G,c,F)

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

lam =      lam_func_2(F,c,d1,d2,d3,G,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot,xg,xg_dot);
ddots = sol_funs_2(F,c,d1,d2,d3,G,lam,m1,m2,m3,m4,theta1,theta2,theta3,theta1_dot,theta2_dot,theta3_dot,xg);

Xdot = [xg_dot, yg_dot, theta1_dot, theta2_dot, theta3_dot, ddots]';
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
