%%Plotting


% thetas = [0 0 0]; %theta's of each pendulum with repsect to vertical
ds = [5 5 5]; % d's of pendulums. Get 3.

% t1 = ones(20,1);
% t2 = zeros(20,1);
% t3 = ones(20,1);
% 
% t1 = gibb(20,0);
% t2 = gibb(20,1);
% t3 = gibb(20,2);
% 
% xg_0 = zeros(20,1);
% yg_0 = zeros(20,1);
% 
% xg_0 = linspace(0,5,20)
% yg_0 = linspace(0,-2,20)

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

for i = 1:length(t1)
    xlim([xg_0(i)-10 xg_0(i)+10])
    ylim([yg_0(i)-20 yg_0(i)+10])
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
    
    pause(0.1)
    
    
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


