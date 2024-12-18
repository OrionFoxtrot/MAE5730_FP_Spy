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

funx0 = linspace(0,60,1e3);
funx1 = linspace(60,100,1e3);
funx2 = linspace(100,140,1e3);
funy0 = (funx0./10-3).^2+1;
funy1 = (-1/10*funx1+16);
funy2 = -(funx2/10-12).^2+10;
plot(funx0,funy0);
plot(funx1,funy1);
plot(funx2,funy2);
xlim([-5 30])
ylim([-10 15])
% xline(0);
% xline(60)
% xline(100)
% xline(140)
xlabel('x')
ylabel('y')
grid on

break_point = -5;
for i = 1:length(t1)
    xlim([xg_0(i)-10 xg_0(i)+10])
    ylim([yg_0(i)-20 yg_0(i)+10])
    
    %CG_0 CG OF SLEDGE
    CG_0 = [xg_0(i), yg_0(i)];

    %CG_1 CG OF LINK 1
    CG_1 = get_next_g(CG_0,ds(1)/2,t1(i));
    link1_end = get_next_g(CG_0,ds(1),t1(i));

    %CG_2 CG OF LINK 2
    CG_2 = get_next_g(link1_end,ds(2)/2,t2(i));
    link2_end = get_next_g(link1_end,ds(1),t2(i));

    %CG_3 CG OF LINK 3
    CG_3 = get_next_g(link2_end,ds(3)/2,t3(i));
    link3_end = get_next_g(link2_end,ds(1),t3(i));

    h1 = plot(CG_0(1),CG_0(2),'o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],MarkerSize=18,DisplayName='CG_0 (Sledge+Head)');
    h2 = plot(CG_1(1),CG_1(2),'o','MarkerFaceColor',[0.4940, 0.1840, 0.5560],MarkerSize=15,DisplayName='CG_1');
    h3 = plot(CG_2(1),CG_2(2),'o','MarkerFaceColor',[0.4660, 0.6740, 0.1880],MarkerSize=15,DisplayName='CG_2');
    h4 = plot(CG_3(1),CG_3(2),'o','MarkerFaceColor',[0.9290, 0.6940, 0.1250],MarkerSize=15,DisplayName='CG_3 (Toes)');

    h5 = plot([CG_0(1),link1_end(1)],[CG_0(2), link1_end(2)],'color','r','LineWidth',2,DisplayName='Link 1');
    h6 = plot([link1_end(1),link2_end(1)],[link1_end(2),link2_end(2)],'color','g','LineWidth',2,DisplayName='Link 2');
    h7 = plot([link2_end(1),link3_end(1)],[link2_end(2),link3_end(2)],'color','b','LineWidth',2, DisplayName='Link 3');
    
    legend([h1 h2 h3 h4 h5 h6 h7],Location='southeast',FontSize=20);
    
    
    title(sprintf("Animation for EL method at index %.2f",i),'FontSize', 30)
    pause(0.1)

    if(break_point == 0)
        text(-3,10,'$$(X_1, Y_1)$$',Interpreter='latex',FontSize=25)
        text(-3,7.5,'$$(X_2, Y_2)$$',Interpreter='latex',FontSize=25)
        text(-3,2.5,'$$(X_3, Y_3)$$',Interpreter='latex',FontSize=25)
        text(-3,-2.5,'$$(X_4, Y_4)$$',Interpreter='latex',FontSize=25)
    end
    if(break_point>=0)
        if(xg_0(i)>break_point)
            break
        end
    end
    

    if i==length(t1)
        break
    end
 
    delete([h1 h2 h3 h4 h5 h6 h7])
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
