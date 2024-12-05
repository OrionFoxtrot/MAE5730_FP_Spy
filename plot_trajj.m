function traj = plot_trajj(X,method)

    if method == 1
        xg_0 = X(:,1);
        yg_0 = X(:,2);
        t1 = X(:,3);
        t2 = X(:,4);
        t3 = X(:,5);
    end
    if method == 2
        t1 = X(:,5);
        t2 = X(:,8);
        t3 = X(:,11);
        xg_0 = X(:,1);
        yg_0 = X(:,2);
    end
    ds = [5 5 5];
    cg1 = [];
    cg2 = [];
    cg3 = [];
    for i = 1:length(t1)
   
        plot(xg_0(i), yg_0(i),DisplayName='G_0 (Head)');
    
        CG_1 = get_next_g([xg_0(i), yg_0(i)],ds(1), t1(i));
        cg1 = [cg1 CG_1];
     
    
        CG_2 = get_next_g(CG_1,ds(2), t2(i));
        cg2 = [cg2 CG_2];
    
        CG_3 = get_next_g(CG_2,ds(3), t3(i));
        cg3 = [cg3 CG_3];
     
    end
    figure
    hold on
    plot(xg_0,yg_0,'MarkerFaceColor','r')
    plot(cg1(1,:), cg1(2,:),'MarkerFaceColor','g')
    plot(cg2(1,:), cg2(2,:),'MarkerFaceColor','b')
    plot(cg3(1,:), cg3(2,:),'MarkerFaceColor','y')
    legend("trolley CG", "Link 1 CG", "Link 2 CG", "Link 3 CG");
    title("trajectory of Trolley");
    
    hold off
  
end

function next_g = get_next_g(Gn,d,theta)
nx = Gn(1) + d*sin(theta);
ny = Gn(2) - d*cos(theta);
next_g = [nx,ny]';
end