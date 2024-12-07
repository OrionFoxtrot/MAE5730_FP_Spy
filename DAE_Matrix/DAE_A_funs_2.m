function A = DAE_A_funs_2(d1,d2,d3,ig1,ig2,ig3,m1,m2,m3,m4,theta1,theta2,theta3,x1)
%DAE_A_funs_2
%    A = DAE_A_funs_2(D1,D2,D3,IG1,IG2,IG3,M1,M2,M3,M4,THETA1,THETA2,THETA3,X1)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Dec-2024 00:37:13

t2 = cos(theta1);
t3 = cos(theta2);
t4 = cos(theta3);
t5 = sin(theta1);
t6 = sin(theta2);
t7 = sin(theta3);
t14 = -m1;
t15 = -m2;
t16 = -m3;
t17 = -m4;
t22 = x1./5.0e+1;
t8 = d1.*t2;
t9 = d2.*t3;
t10 = d3.*t4;
t11 = d1.*t5;
t12 = d2.*t6;
t13 = d3.*t7;
t18 = t8.*2.0;
t19 = t11.*2.0;
t20 = -t11;
t21 = -t12;
mt1 = [t14,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,t22-1.2e+1./5.0,0.0,0.0,t14,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,t15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-ig1,0.0,0.0,0.0,0.0,0.0,0.0,t8,t11,t18,t19,t18,t19,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-ig2,0.0,0.0,0.0,0.0,0.0,t9,t12,t9.*2.0,t12.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t17,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
mt2 = [0.0,t17,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-ig3,0.0,0.0,0.0,0.0,t10,t13,0.0,0.0,1.0,0.0,-1.0,0.0,t8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,t8,-1.0,0.0,t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t20,0.0,1.0,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,t9,-1.0,0.0,t10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,t21,0.0,1.0,-t13,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
mt3 = [-t22+1.2e+1./5.0];
A = reshape([mt1,mt2,mt3],19,19);
end
