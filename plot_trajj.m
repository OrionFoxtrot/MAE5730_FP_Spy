function traj = plot_trajj(X,method)

    if method == 1
        xg_0 = X(:,1);
        yg_0 = X(:,2);
        t1 = X(:,3);
        t2 = X(:,4);
        t3 = X(:,5);
        met = "Euler-Lagrange";
    end
    if method == 2
        t1 = X(:,5);
        t2 = X(:,8);
        t3 = X(:,11);
        xg_0 = X(:,1);
        yg_0 = X(:,2);
        met = "DAE";
    end
    ds = [5 5 5];
    cg1 = [];
    cg2 = [];
    cg3 = [];
    for i = 1:length(t1)
  
        CG_1 = get_next_g([xg_0(i), yg_0(i)],ds(1), t1(i));
        cg1 = [cg1 CG_1];
    
        CG_2 = get_next_g(CG_1,ds(2), t2(i));
        cg2 = [cg2 CG_2];
    
        CG_3 = get_next_g(CG_2,ds(3), t3(i));
        cg3 = [cg3 CG_3];
     
    end
 
    hold on
    funx0 = linspace(0,60,1e3);
    funx1 = linspace(60,100,1e3);
    funx2 = linspace(100,140,1e3);
    funy0 = (funx0./10-3).^2+1;
    funy1 = (-1/10*funx1+16);
    funy2 = -(funx2/10-12).^2+10;
    %plot(funx0,funy0,funx1,funy1,funx2,funy2,DisplayName="Course")

    plot(xg_0,yg_0,'r',DisplayName="Trolley CG")
    plot(cg1(1,:), cg1(2,:),'g',DisplayName="Link 1 CG")
    plot(cg2(1,:), cg2(2,:),'b',DisplayName="Link 2 CG")
    plot(cg3(1,:), cg3(2,:),'k',DisplayName="Link 3 CG")
    legend()
    title(sprintf("Trajectory of trolley with %s method",met));
    
    hold off
    f = gcf;
    exportgraphics(f,sprintf('Trajectories/%s.png',met),'Resolution',300)
  
end

function next_g = get_next_g(Gn,d,theta)
nx = Gn(1) + d*sin(theta);
ny = Gn(2) - d*cos(theta);
next_g = [nx,ny]';
end