%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift and drag coefficients with Alpha——随攻角变化的升阻力系数
clear all; clc;
% alpha=atan2(-w_y,w_z);    %Aerodynamic angle of attack  %这个要求事先求出绕各轴的翅运动角速率
beta=pi/4; f=0.25; %\Hz
delta=pi/4;    %相差delta:  提前模式:delta=pi/4; 延迟模式:delta=-pi/4; 对称模式:delta=0;
t=linspace(0,10,200);
alpha=(beta*sin(2*pi*f*t+delta))*180/pi;
% figure(1)
% plot(t,alpha,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 先使用函数fzero,求零点值，再使用函数find找到其位置索引
% syms  t
% beta=pi/4; f=0.25; %\Hz
% delta=pi/4;    %相差delta:  提前模式，delta=-pi/4为延迟模式; delta=0，对称模式;
% alpha=(beta*sin(2*pi*f*t+delta))*180/pi
% % Result
% % alpha = 45*sin(pi/4 + (pi*t)/2)
options=optimset('Display','off');
t1_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[3,4],options); 
t2_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[5,6],options);
% Result:   t1_0=3.5000;   t2_0=5.5000;
% indxx1=find(t==1.5)  %取等号值'=='
indx1=find(t>3.45 & t<=3.5);  %条件语句
indx2=find(t>5.45 & t<=5.5);
% Result:   indx1=70;   indx2=110;
alpha_2=alpha(1,indx1:indx2);      % 半个周期    %  alpha_2∈(0°,45°)
% figure(2)
% plot(t(1,indx1:indx2),alpha_2,'ro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 求出alpha_2的最大值和最小值及其位置索引
[alpha_2_max, locat_max]=max(alpha_2);   % 输出:  alpha_2_max =44.9716;  locat_max =22;
% [alpha_2_min locat_min]=min(alpha_2);  % 输出:  alpha_2_min =-2.3078;    locat_min =1;
m=length(alpha_2);                                                     %m =41
alpha_3=alpha_2(1,locat_max:m)+45;     % size(alpha_3)  % (1*20)                  %  alpha_3∈(45°,90°)
[alpha_new,index]=sort(alpha_3);           %按升序排列，返回两个数据
%这里将第一项alpha_2(1,1)设为0;三组数据重新构成新的向量——其实可以直接通过linspace函数生成(0°,90°)的数据序列
alpha_4=[0, alpha_2(1,2:locat_max-1), alpha_new];  % size(alpha_4)% (1*41)   %  alpha_4∈(0°,90°)
% figure(3)
% plot(t(1,indx1:indx2),alpha_4,'m-d')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 各种升阻力系数的对比分析
%% (1) 1999-Science-MH Dickinson——三维实验测试升阻力系数
% C_rot_theo=pi*(0.75-x_0nd(r));　  % 转动环量气动力系数
C_L_dickinson=0.225+1.58*sin((2.13*alpha_4-7.2)*pi/180);    % 角度alpha以度数表示,但是三角函数是针对弧度制的哦
C_D_dickinson=1.92-1.55*cos((2.04*alpha_4-9.82)*pi/180);   % 角度alpha以度数表示,但是三角函数是针对弧度制的哦
% figure(4)
% plot(alpha_4,C_L_dickinson,'r-',alpha_4,C_D_dickinson,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 2005-ACC-Khan
% C_r=0.75-x_0nd(r);    % C_r=0.75-d_0(r);　　% 转动环量气动力系数
C_t_khan=7*abs(alpha_4*pi/180)/pi;  % 需要转换成升阻力系数
C_L_khan=C_t_khan.*sin(2*alpha_4*pi/180); 
C_D_khan=C_L_khan.*tan(alpha_4*pi/180); 
% figure(5)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_khan,'r-',alpha_4,C_D_khan,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 2004-JEB-Wang ZJ & 2005-JFM-Andersen_b & 2007-JFM-Bergou
% C_R=pi;    % 2005-JFM-Andersen_b
% 文献2007-JFM-Bergou中，针对蜻蜓前后翅，对比准稳态和CFD结果的差值，给出了分别是
% C_R_forewing=2.58; C_R_hindwing=0.9; 
% C_T=1.2;                 % 2005-JFM-Andersen_b
% A=1.4;   B=1.0;      % 2005-JFM-Andersen_b
% C_L_Andersen=C_T*sin(2*alpha_4*pi/180); 
% C_D_Andersen=A-B*cos(2*alpha_4*pi/180);
A=1.2; B=1.4;   C=1.0;
C_L_wang=A*sin(2*alpha_4*pi/180);       % 由二维CFD计算数据拟合而得,
C_D_wang=B-C*cos(2*alpha_4*pi/180);  % 由二维CFD计算数据拟合而得,可能是针对震颤跌宕至稳态下落的卡片
% figure(6)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_wang,'r-',alpha_4,C_D_wang,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) 2005-JFM-Andersen_a
% C_R=1.0;               % C_R∈(1.0,1.4)  % 2005-JFM-Andersen_a
% C_T=1.0;                   % 2005-JFM-Andersen_a
% C_d_alpha0=0.13;              % 2005-JFM-Andersen_a  % C_d0∈(0.1,0.3)
% C_d_alphapi_2=2.0;               % 2005-JFM-Andersen_a
C_T=1.833;                  % 2007-JFM-Berman  % Wang ZJ拟合1999-Science-Dickinson数据
C_d_alpha0=0.21;                  % 2007-JFM-Berman % Wang ZJ拟合1999-Science-Dickinson数据
C_d_alphapi_2=3.35;                 % 2007-JFM-Berman % Wang ZJ拟合1999-Science-Dickinson数据  
C_L_andersen=C_T*sin(2*alpha_4*pi/180); 
C_D_andersen=C_d_alpha0*(cos(alpha_4*pi/180)).^2+C_d_alphapi_2*(sin(alpha_4*pi/180)).^2;  % 针对震颤跌宕的下落纸板
% C_L_andersen=C_D_andersen./tan(alpha_4*pi/180);  % 升力系数;参考2014-JRSI-Nabawy
% figure(7)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_andersen,'r-',alpha_4,C_D_andersen,'b-') 
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) 2010-JFM-Whitney
% Here mentioned coefficients is for the situation of the Reynolds number (≈100), 
% geometry and flapping kinematics of a typical Drosophila wing, Dickinson et al. (1999)
% C_R=1.55;   % 转动环量气动力系数随着转动速率和转动轴的位置变化 
C_Lmax=1.8;   C_Dmax=3.4;     C_D0=0.4;  
C_L_whitney=C_Lmax*sin(2*alpha_4*pi/180);  
C_D_whitney=(C_Dmax+C_D0)/2-(C_Dmax-C_D0)/2*cos(2*alpha_4*pi/180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) 2006-IEEE TR-Deng Xinyan & 2010-EAE-RJ Wood
% C_rot_theo=2*pi*(0.75-x_0nd(r));　  % 转动环量气动力系数
% C_N_deng=3.4*sin(abs(alpha_4*pi/180));     % 2010-EAE-RJ Wood           % 来自1999-Science-Dickinson的实验数据
% C_T_deng=0.4*(cos(2*alpha_4*pi/180)).^2.*(alpha_4>-45 & alpha_4<45);   % 当alpha∈(-pi/4,pi/4), or else C_T=0;
C_N_deng=3.4*sin(alpha_4*pi/180);         % 2006-IEEE TR-Deng Xinyan % 来自1999-Science-Dickinson的实验数据
C_T_deng=0.4*(cos(2*alpha_4*pi/180)).^2.*(alpha_4>=0 & alpha_4<=45);   % 当alpha∈[0,pi/4], or else C_T=0;
C_L_deng=C_N_deng.*cos(alpha_4*pi/180)-C_T_deng.*sin(alpha_4*pi/180);
C_D_deng=C_N_deng.*sin(alpha_4*pi/180)+C_T_deng.*cos(alpha_4*pi/180);
% figure(9)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_N_deng,'k-',alpha_4,C_T_deng,'r-',alpha_4,C_L_deng,'g-',alpha_4,C_D_deng,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) 2010-DSCC-Katie_Byl
% 简化的目的是：
% These simplified curves also ensure both smoothness and symmetry (as theoretically required) near alpha_4=0;
% aid the reader in understanding the basic, theoretical relationships between the passive dynamics and output control forces more intuitively.
C_L_byl=1.8*sin(2*alpha_4*pi/180);
C_D_byl=1.8*(1-cos(2*alpha_4*pi/180));
% C_T_byl=0;  
C_N_byl=3.6*sin(alpha_4*pi/180);   % C_N_byl=sqrt(C_L_byl.^2+C_D_byl.^2);
% figure(10)
% plot(alpha_4,C_L_byl,'k-',alpha_4,C_D_byl,'b-',alpha_4,C_N_byl,'b-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) 波尔豪森采用前缘吸力比拟法建立的升力系数计算模型
% C_L_polhamus=K_p*sin(alpha)*(cos(alpha)).^2+K_v*(cos(alpha))*sin(alpha).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) 2014-JRSI-Nabawy
C_La2d=2*pi;    % 单位是rad^(-1), 这里2*pi (rad^(-1))= 0.11 (deg^(-1));   
% C_La2d=5.15662;  % 单位是rad^(-1), 针对典型昆虫雷诺数的平板, 建议使用 0.09 (deg^(-1))
% The parameter E is the edge correction proposed by Jones [42] for the lifting line theory 
% and is evaluated as the quotient of the wing semi-perimeter to its length[42,43].
% lambda=0.627; % 人为引入一个缩放因子，以便于实验测试而得的气动力系数进行对比
lambda=0.755;
E=lambda*1.146071666;   % E=C/2/R;
% The parameter k is the so-called 'k-factor' included to correct for the difference in efficiency between the assumed ideal
% uniform downwash distribution and real downwash distribution [35,40,44]. For this work, the k-factor required within
% equation (2.11) will be estimated using the induced power factor expression of hovering actuator disc models [45].
k=1.51;  % 针对果蝇
% for the span of a single  wing 
AR=3.40158;     % 注意这里的展弦比AR∈(3,5)是比较好的。% 针对果蝇
C_La_Ny=C_La2d/(E+k*C_La2d/(pi*AR));     % 由普朗特升力线理论获得的三维翅膀升力曲线斜率  % 输出: C_La_Ny =2.7506;
% C_L_nabawy=(0.5*C_La2d/(E+k*C_La2d/(pi*AR)))*sin(2*alpha_4);
C_L_nabawy=0.5*C_La_Ny*sin(2*alpha_4*pi/180);   % 与2009-JEB-Rotational accelerations 的Fig.7图实测数据对比符合很好
% C_D_nabawy=C_Lnabawy.*tan(alpha_4);   % C_Dnabawy=2*C_T*(sin(alpha))^2=C_La_Ny*(sin(alpha))^2;  
C_D_nabawy=C_La_Ny*(sin(alpha_4*pi/180)).^2;  
% 下面是采用前缘吸力比拟法建立的波尔豪森升力系数模型
% C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha)).*(cos(alpha)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha));
C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha_4*pi/180)); 
% figure(11)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_nabawy,'r-',alpha_4,C_D_nabawy,'b-',alpha_4,C_L_polhamus1,'g-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) 2014-AST-Taha
% a0 is the lift curve slope of the two-dimensional airfoil section, e.g. it is equal to 2π for a flat plate or a very thin cambered shape.
% For conventional airfoils,it could be determined from lift curves such as the ones presented by Abbott and Doenhoff[1].
a_0=2*pi;
%  AR being based on one wing: AR=R^2/S;
AR=2.8;  % 注意这里的展弦比AR＜3是比较好的。 
C_La_Ta=pi*AR/(1+sqrt((pi*AR/a_0)^2+1));     % 由扩展升力线获得的三维翅膀升力曲线斜率   % 输出:C_La_Ta =3.2334;
% C_L_taha=(0.5*pi*AR/(1+sqrt((pi*AR/a_0)^2+1))).*sin(2*alpha_4*pi/180);
C_L_taha=0.5*C_La_Ta*sin(2*alpha_4*pi/180);
C_D_taha=C_L_taha.*tan(alpha_4*pi/180);
% 下面是采用前缘吸力比拟法建立的波尔豪森升力系数模型
% % C_L_polhamus2=C_La_Ta*sin(alpha).*(cos(alpha)).^2+(C_La_Ta-C_La_Ta^2/(pi*AR)).*cos(alpha)*(sin(alpha)).^2;  
% C_L_polhamus2=C_La_Ta*sin(alpha_4*pi/180).*(cos(alpha_4*pi/180)).^2+(C_La_Ta-C_La_Ta^2/(pi*AR))*cos(alpha_4*pi/180).*(sin(alpha_4*pi/180)).^2; 
C_L_polhamus2=(0.5*C_La_Ta*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-C_La_Ta/(pi*AR)).*sin(alpha_4*pi/180));
cftool  % sum of sine 谐波函数拟合 or Fourier 级数拟合 平动气动力系数
figure(12)
plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_taha,'r-',alpha_4,C_D_taha,'b',...
       alpha_4,C_L_nabawy,'k-',alpha_4,C_D_nabawy,'k-',alpha_4,C_L_polhamus1,'c-',alpha_4,C_L_polhamus2,'m-')
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha ) and \it C_D (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{D,dickinson}','\itC_{L,taha}','\itC_{D,taha}','\itC_{L,nabawy}',...
            '\itC_{D,nabawy}','\itC_{L,polhamus1}','\itC_{L,polhamus2}')       
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一部分——升力系数的对比分析
figure(13)    %  subplot(224)
C_L=plot(alpha_4,C_L_dickinson,'k-',alpha_4,C_L_khan,'r--',...
                alpha_4,C_L_wang,'r-.',alpha_4,C_L_andersen,'r:',alpha_4,C_L_whitney,'b--',...
                alpha_4,C_L_deng,'b-.',alpha_4,C_L_byl,'b:',alpha_4,C_L_nabawy,'g--',alpha_4,C_L_polhamus1,'g-.',...
                alpha_4,C_L_taha,'g:',alpha_4,C_L_polhamus2,'m--');
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{L,khan}','\itC_{L,wang}','\itC_{L,andersen}','\itC_{L,whitney}',...
             '\itC_{L,deng}','\itC_{L,byl}','\itC_{L,nabawy}','\itC_{L,polhamus1}','\itC_{L,taha}','\itC_{L,polhamus2}')
set(C_L,'LineWidth',2) 
grid on
%axis([xmin,xmax,ymin,ymax])
%% 第二部分——阻力系数的对比分析
figure(14)
C_D=plot(alpha_4,C_D_dickinson,'k-',alpha_4,C_D_khan,'r--',...
                alpha_4,C_D_wang,'r-.',alpha_4,C_D_andersen,'r:',alpha_4,C_D_whitney,'b--',...
                alpha_4,C_D_deng,'b-.',alpha_4,C_D_byl,'b:',alpha_4,C_D_nabawy,'g--',alpha_4,C_D_taha,'g-.');
xlabel('\it\alpha (deg.)')
ylabel('\it C_D (\alpha )')
title('Aerodynamic coefficients of drag \itvs. \alpha \rm for flapping wing')
legend('\itC_{D,dickinson}','\itC_{D,khan}','\itC_{D,wang}','\itC_{D,andersen}','\itC_{D,whitney}',...
              '\itC_{D,deng}','\itC_{D,byl}','\itC_{D,nabawy}','\itC_{D,taha}')
set(C_D,'LineWidth',2)     %参见——《MATLAB原理与工程应用》(第二版)高会声等/译 P182-183
grid on
%% 坐标轴显示区间的设定
v_axis=axis;     %axis([xmin,xmax,ymin,ymax])
%Result: v_axis =  0   100     0     4
% v_axis(1)=0;          %指定x轴的最小值
% v_axis(2)=100;      %指定x轴的最大值
v_axis(3)=-0.5;          %指定y轴的最小值
v_axis(4)=3.5;           %指定y轴的最大值
axis(v_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求出alpha的最大值和最小值及其位置索引
% [alpha_max, locat_max]=max(alpha);
% [alpha_min, locat_min]=min(alpha);
% %Result:
% % alpha_max =134.9996
% % locat_max = 11
% % alpha_min =45.0088
% % locat_min =51
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % alpha=(beta*sin(2*pi*f*t+delta))*180/pi                      %几何攻角的alpha变化规律
%%% Result: alpha = 45*sin(pi/4 + (pi*t)/2)
%%%%%%%%%%%%%%%%%%%%
%求出alpha的最大值和最小值及其位置索引
% options=optimset('Display','off');
% [tmax fmax]=fminbnd(inline('-45*sin(pi/4 + (pi*t)/2)','t'),0,1,options); 
% disp(['Maximum value of alpha in the interval 0<=t<=1 is '  num2str(-fmax)])
% disp(['which occurs at t= '  num2str(tmax)])
% [tmin fmin]=fminbnd(inline('45*sin(pi/4 + (pi*t)/2)','t'),2,3,options); 
% disp(['Minimum value of alpha in the interval 2<=t<=3 is '  num2str(fmin)])
% disp(['which occurs at t= '  num2str(tmin)])
% % Result
% % Maximum value of alpha in the interval 0<=t<=1 is 45
% % which occurs at t= 0.5
% % Minimum value of alpha in the interval 2<=t<=3 is -45
% % which occurs at t= 2.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha_0=pi/2; 
% alpha=(alpha_0+beta*sin(2*pi*f*t+delta))*180/pi
%%% Result: alpha =(90*pi + 45*pi*sin(pi/4 + (pi*t)/2))/pi
%%%%%%%%%%%%%%%%%%%%
%求出alpha的最大值和最小值及其位置索引
% options=optimset('Display','off');
% [tmax fmax]=fminbnd(inline('-(90*pi + 45*pi*sin(pi/4 + (pi*t)/2))/pi','t'),0,1,options); 
% disp(['Maximum value of alpha in the interval 0<=t<=1 is '  num2str(-fmax)])
% disp(['which occurs at t= '  num2str(tmax)])
% [tmin fmin]=fminbnd(inline('(90*pi + 45*pi*sin(pi/4 + (pi*t)/2))/pi','t'),2,3,options); 
% disp(['Minimum value of alpha in the interval 2<=t<=3 is '  num2str(fmin)])
% disp(['which occurs at t= '  num2str(tmin)])
% % Result
% % Maximum value of alpha in the interval 0<=t<=1 is 135
% % which occurs at t= 0.5
% % Minimum value of alpha in the interval 2<=t<=3 is 45
% % which occurs at t= 2.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


