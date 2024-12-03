function traj = plot_traj(X)

    xg_0 = X(:,1);
    yg_0 = X(:,2);
    
    xg_1 = X(:,3);
    yg_1 = X(:,4);

    xg_2 = X(:,6);
    yg_2 = X(:,7);

    xg_3 = X(:,9);
    yg_3 = X(:,10);

    hold on
    plot(xg_0,yg_0);
    plot(xg_1,yg_1);
    plot(xg_2,yg_2);
    plot(xg_3,yg_3);
    legend("trolley CG", "Link 1 CG", "Link 2 CG", "Link 3 CG");
    title("trajectory of Trolley");
    hold off
end

function next_g = get_next_g(Gn,d,theta)
nx = Gn(1) + d*sin(theta);
ny = Gn(2) - d*cos(theta);
next_g = [nx,ny]';
end