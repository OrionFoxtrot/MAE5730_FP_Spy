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

    h1 = plot(CG_0(1),CG_0(2),'o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],MarkerSize=10,DisplayName='CG_0 (Sledge+Head)');
    h2 = plot(CG_1(1),CG_1(2),'o','MarkerFaceColor',[0.4940, 0.1840, 0.5560],MarkerSize=10,DisplayName='CG_1');
    h3 = plot(CG_2(1),CG_2(2),'o','MarkerFaceColor',[0.4660, 0.6740, 0.1880],MarkerSize=10,DisplayName='CG_2');
    h4 = plot(CG_3(1),CG_3(2),'o','MarkerFaceColor',[0.9290, 0.6940, 0.1250],MarkerSize=10,DisplayName='CG_3 (Toes)');

    h5 = plot([CG_0(1),link1_end(1)],[CG_0(2), link1_end(2)],'color','r','LineWidth',1,DisplayName='Link 1');
    h6 = plot([link1_end(1),link2_end(1)],[link1_end(2),link2_end(2)],'color','g','LineWidth',1,DisplayName='Link 2');
    h7 = plot([link2_end(1),link3_end(1)],[link2_end(2),link3_end(2)],'color','b','LineWidth',1, DisplayName='Link 3');
    
    legend([h1 h2 h3 h4 h5 h6 h7],Location='southeast',FontSize=20);

    title(sprintf("Animation for DAE method at index %.2f",i),'FontSize', 24)
    
    pause(0.01)
    
    % if(xg_0(i)>20)
    %     break
    % end
    
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


