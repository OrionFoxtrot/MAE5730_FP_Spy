%%Plotting

%using X0's unpack
xg_0 = X(:,1);
yg_0 = X(:,2);

xg_1 = X(:,3);
yg_1 = X(:,4);

xg_2 = X(:,6);
yg_2 = X(:,7);

xg_3 = X(:,9);
yg_3 = X(:,10);

% xg_0 = [X(:,1); X1(:,1)];
% yg_0 = [X(:,2); X1(:,2)];
% 
% xg_1 = [X(:,3); X1(:,3)];
% yg_1 = [X(:,4); X1(:,4)];
% 
% xg_2 = [X(:,6); X1(:,6)];
% yg_2 = [X(:,7); X1(:,7)];
% 
% xg_3 = [X(:,9); X1(:,9)];
% yg_3 = [X(:,10); X1(:,10)];

%overlay stuff
l = length(xg_0);
funx = linspace(min(xg_0),max(xg_0));
funy = 2*cos(funx/3);

close all
hold on
plot(funx,funy);
axis equal
xlim([-10 10])
ylim([-20 10])

for i = 1:length(xg_0)
    xlim([xg_0(i)-10 xg_0(i)+10])
    ylim([yg_0(i)-20 yg_0(i)+10])

    h1 = plot(xg_0(i), yg_0(i), 'o','MarkerFaceColor','r',MarkerSize=15, DisplayName='G_0 (Head)');
    h2 = plot(xg_1(i), yg_1(i), 'o','MarkerFaceColor','b',MarkerSize=10,DisplayName='G_1');
    h3 = plot(xg_2(i), yg_2(i), 'o','MarkerFaceColor','g',MarkerSize=10,DisplayName='G_2');
    h4 = plot(xg_3(i), yg_3(i), 'o','MarkerFaceColor','y',MarkerSize=10,DisplayName='G_3 (Toes)');

    h5 = plot( [xg_0(i),xg_1(i)], [yg_0(i),yg_1(i)],'MarkerFaceColor','k');
    h6 = plot( [xg_1(i),xg_2(i)], [yg_1(i),yg_2(i)],'MarkerFaceColor','k');
    h7 = plot( [xg_2(i),xg_3(i)], [yg_2(i),yg_3(i)],'MarkerFaceColor','k');

    

    legend([h1,h2,h3,h4])
    title(sprintf("Index %.2f",i))
    
    pause(0.05)
    
    
    if i==length(xg_0)
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


