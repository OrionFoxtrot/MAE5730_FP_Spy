%%Plotting


% thetas = [0 0 0]; %theta's of each pendulum with repsect to vertical
ds = [5 5 5]; % d's of pendulums. Get 3.

%using X0's unpack
xg_0 = X(:,1);
yg_0 = X(:,2);
t1 = X(:,3);
t2 = X(:,4);
t3 = X(:,5);

close all
hold on
axis equal
xlim([-10 10])
ylim([-20 10])

% funx = linspace(min(xg_0)-5,max(xg_0)+5,1e5);
%funy = funx; %y=x
%funy = (funx-3).^2+1; %y = (x-3)^2+1
% funy = (funx/10-3).^2+1; %y = (x/10-3)^2+1
%funy = cos(funx);
%xline(2*pi)
% plot(funx,funy)

funx0 = linspace(0,60,1e3);
funx1 = linspace(60,100,1e3);
funx2 = linspace(100,140,1e3);
funy0 = (funx0./10-3).^2+1;
funy1 = (-1/10*funx1+16);
funy2 = -(funx2/10-12).^2+10;
plot(funx0,funy0);
plot(funx1,funy1);
plot(funx2,funy2);
xlim([-10 30])
ylim([-15 15])
xline(0);
xline(60)
xline(100)
xline(140)
xlabel('x')
ylabel('y')

for i = 1:length(t1)
    % xlim([xg_0(i)-10 xg_0(i)+10])
    % ylim([yg_0(i)-20 yg_0(i)+10])
    h1 = plot(xg_0(i), yg_0(i), 'o','MarkerFaceColor','r',MarkerSize=15, DisplayName='G_0 (Head)');

    CG_1 = get_next_g([xg_0(i), yg_0(i)],ds(1), t1(i));
    h2 = plot(CG_1(1), CG_1(2), 'o','MarkerFaceColor','b',MarkerSize=10,DisplayName='G_1');

    CG_2 = get_next_g(CG_1,ds(2), t2(i));
    h3 = plot(CG_2(1), CG_2(2), 'o','MarkerFaceColor','g',MarkerSize=10,DisplayName='G_2');

    CG_3 = get_next_g(CG_2,ds(3), t3(i));
    h4 = plot(CG_3(1), CG_3(2), 'o','MarkerFaceColor','y',MarkerSize=10,DisplayName='G_3 (Toes)');

    h5 = plot( [xg_0(i),CG_1(1)], [yg_0(i),CG_1(2)],'MarkerFaceColor','k');
    h6 = plot( [CG_1(1),CG_2(1)], [CG_1(2),CG_2(2)],'MarkerFaceColor','k');
    h7 = plot( [CG_2(1),CG_3(1)], [CG_2(2),CG_3(2)],'MarkerFaceColor','k');

    legend([h1,h2,h3,h4])
    
    title(sprintf("Animation for EL method at index %.2f",i))
    pause(0.05)
    
    if(xg_0(i)>0)
        break
    end

    if i==length(t1)
        break
    end
    delete([h5 h6 h7])
    delete([h1 h2 h3 h4])
end

function next_g = get_next_g(Gn,d,theta)
nx = Gn(1) + d*sin(theta);
ny = Gn(2) - d*cos(theta);
next_g = [nx,ny]';
end

function gib_sin = gibb(n,off)
z = zeros(1,n);
for i = 1:n
    z(i) = sin(i+off);
end
gib_sin = z;
end


%% Energy
% X0_0 = [0, 10, 0,0,0];
% X0_1 = [0,0,0,0,0]; %xg_dot yg_dot, theta1_dot, theta2_dot, theta3_dot
% X0 = [X0_0, X0_1];
xg = X(:,1);
yg = X(:,2);
theta1 = X(:,3);
theta2 = X(:,4);
theta3 = X(:,5);

xg_dot = X(:,6);
yg_dot = X(:,7);
theta1_dot = X(:,8);
theta2_dot = X(:,9);
theta3_dot = X(:,10);

m1 = 30;
m2 = 15;
m3 = 5;
m4 = 2;

d1 = 5;
d2 = 5;
d3 = 5;

g = 9.81;

Ig1 = 1/12*m2*(2*d1).^2;
Ig2 = 1/12*m3*(2*d2).^2;
Ig3 = 1/12*m4*(2*d3).^2;


vg = sqrt(xg_dot.^2+yg_dot.^2);
vg2 = sqrt((xg_dot + d1.*cos(theta1).*theta1_dot).^2 + (yg_dot + d1.*sin(theta1).*theta1_dot).^2);
vg3 = sqrt((xg_dot + 2.*d1.*cos(theta1).*theta1_dot + d2.*cos(theta2).*theta2_dot).^2 + (yg_dot + 2.*d1.*sin(theta1).*theta1_dot + d2.*sin(theta2).*theta2_dot).^2);
vg4 = sqrt((xg_dot + 2.*d1.*cos(theta1).*theta1_dot + 2.*d2.*cos(theta2).*theta2_dot + d3.*cos(theta3).*theta3_dot).^2 + (yg_dot + 2.*d1.*sin(theta1).*theta1_dot + 2.*d2.*sin(theta2).*theta2_dot + d3.*sin(theta3).*theta3_dot).^2);

KE = 1/2.* ( m1.*vg.^2 ...
    +m2.*vg2.^2+Ig1.*theta1_dot.^2 ...
    +m3.*vg3.^2+Ig2.*theta2_dot.^2 ...
    +m4.*vg4.^2+Ig3.*theta3_dot.^2 );
PE = g.* (m1.*yg ...
    + m2.*(yg-d1.*cos(theta1))...
    + m3.*(yg-2.*d1.*cos(theta1)-d2.*cos(theta2))...
    + m4.*(yg -2.*d1.*cos(theta1)-2.*d2.*cos(theta2)-d3.*cos(theta3)));
% KE = KE-min(KE);
% PE = PE-min(PE);
close all
hold on
plot(KE,'r','LineWidth',2)
plot(PE,'g','LineWidth',2)
plot(KE+PE,'--','color','b','LineWidth',2)
legend('Kinetic Energy', 'Potential Energy', 'Total Energy','FontSize', 24)
title('Energy in system','FontSize', 24)
xlabel('time (s)','FontSize', 24)
ylabel('Energy (J)','FontSize', 24)
