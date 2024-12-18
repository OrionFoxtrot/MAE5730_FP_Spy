%%Plotting

%using X0's unpack
ds = [5,5,5];
t1 = X(:,5);
t2 = X(:,8);
t3 = X(:,11);
xg_0 = X(:,1);
yg_0 = X(:,2);

close all
hold on

%overlay stuff
funx0 = linspace(0,60,1e3);
funx1 = linspace(60,100,1e3);
funx2 = linspace(100,140,1e3);
funy0 = (funx0./10-3).^2+1;
funy1 = (-1/10*funx1+16);
funy2 = -(funx2/10-12).^2+10;
plot(funx0,funy0);
plot(funx1,funy1);
plot(funx2,funy2);

axis equal
xlim([-10 60])
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
    h1 = plot(xg_0(i), yg_0(i), 'o','MarkerFaceColor','r',MarkerSize=30, DisplayName='G_0 (Head)');

    CG_1 = get_next_g([xg_0(i), yg_0(i)],ds(1), t1(i));
    h2 = plot(CG_1(1), CG_1(2), 'o','MarkerFaceColor','b',MarkerSize=20,DisplayName='G_1');

    CG_2 = get_next_g(CG_1,ds(2), t2(i));
    h3 = plot(CG_2(1), CG_2(2), 'o','MarkerFaceColor','g',MarkerSize=20,DisplayName='G_2');

    CG_3 = get_next_g(CG_2,ds(3), t3(i));
    h4 = plot(CG_3(1), CG_3(2), 'o','MarkerFaceColor','y',MarkerSize=20,DisplayName='G_3 (Toes)');

    h5 = plot( [xg_0(i),CG_1(1)], [yg_0(i),CG_1(2)],'MarkerFaceColor','k');
    h6 = plot( [CG_1(1),CG_2(1)], [CG_1(2),CG_2(2)],'MarkerFaceColor','k');
    h7 = plot( [CG_2(1),CG_3(1)], [CG_2(2),CG_3(2)],'MarkerFaceColor','k');

    legend([h1,h2,h3,h4],'FontSize', 24,'Location','southeast')
    title(sprintf("Animation for DAE method at index %.2f",i),'FontSize', 24)
    
    pause(0.01)
    
    if(xg_0(i)>20)
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


