%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift and drag coefficients with Alpha�����湥�Ǳ仯��������ϵ��
clear all; clc;
% alpha=atan2(-w_y,w_z);    %Aerodynamic angle of attack  %���Ҫ����������Ƹ���ĳ��˶�������
beta=pi/4; f=0.25; %\Hz
delta=pi/4;    %���delta:  ��ǰģʽ:delta=pi/4; �ӳ�ģʽ:delta=-pi/4; �Գ�ģʽ:delta=0;
t=linspace(0,10,200);
alpha=(beta*sin(2*pi*f*t+delta))*180/pi;
% figure(1)
% plot(t,alpha,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��ʹ�ú���fzero,�����ֵ����ʹ�ú���find�ҵ���λ������
% syms  t
% beta=pi/4; f=0.25; %\Hz
% delta=pi/4;    %���delta:  ��ǰģʽ��delta=-pi/4Ϊ�ӳ�ģʽ; delta=0���Գ�ģʽ;
% alpha=(beta*sin(2*pi*f*t+delta))*180/pi
% % Result
% % alpha = 45*sin(pi/4 + (pi*t)/2)
options=optimset('Display','off');
t1_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[3,4],options); 
t2_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[5,6],options);
% Result:   t1_0=3.5000;   t2_0=5.5000;
% indxx1=find(t==1.5)  %ȡ�Ⱥ�ֵ'=='
indx1=find(t>3.45 & t<=3.5);  %�������
indx2=find(t>5.45 & t<=5.5);
% Result:   indx1=70;   indx2=110;
alpha_2=alpha(1,indx1:indx2);      % �������    %  alpha_2��(0��,45��)
% figure(2)
% plot(t(1,indx1:indx2),alpha_2,'ro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���alpha_2�����ֵ����Сֵ����λ������
[alpha_2_max, locat_max]=max(alpha_2);   % ���:  alpha_2_max =44.9716;  locat_max =22;
% [alpha_2_min locat_min]=min(alpha_2);  % ���:  alpha_2_min =-2.3078;    locat_min =1;
m=length(alpha_2);                                                     %m =41
alpha_3=alpha_2(1,locat_max:m)+45;     % size(alpha_3)  % (1*20)                  %  alpha_3��(45��,90��)
[alpha_new,index]=sort(alpha_3);           %���������У�������������
%���ｫ��һ��alpha_2(1,1)��Ϊ0;�����������¹����µ�����������ʵ����ֱ��ͨ��linspace��������(0��,90��)����������
alpha_4=[0, alpha_2(1,2:locat_max-1), alpha_new];  % size(alpha_4)% (1*41)   %  alpha_4��(0��,90��)
% figure(3)
% plot(t(1,indx1:indx2),alpha_4,'m-d')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����������ϵ���ĶԱȷ���
%% (1) 1999-Science-MH Dickinson������άʵ�����������ϵ��
% C_rot_theo=pi*(0.75-x_0nd(r));��  % ת������������ϵ��
C_L_dickinson=0.225+1.58*sin((2.13*alpha_4-7.2)*pi/180);    % �Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�Ŷ
C_D_dickinson=1.92-1.55*cos((2.04*alpha_4-9.82)*pi/180);   % �Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�Ŷ
% figure(4)
% plot(alpha_4,C_L_dickinson,'r-',alpha_4,C_D_dickinson,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 2005-ACC-Khan
% C_r=0.75-x_0nd(r);    % C_r=0.75-d_0(r);����% ת������������ϵ��
C_t_khan=7*abs(alpha_4*pi/180)/pi;  % ��Ҫת����������ϵ��
C_L_khan=C_t_khan.*sin(2*alpha_4*pi/180); 
C_D_khan=C_L_khan.*tan(alpha_4*pi/180); 
% figure(5)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_khan,'r-',alpha_4,C_D_khan,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 2004-JEB-Wang ZJ & 2005-JFM-Andersen_b & 2007-JFM-Bergou
% C_R=pi;    % 2005-JFM-Andersen_b
% ����2007-JFM-Bergou�У��������ǰ��ᣬ�Ա�׼��̬��CFD����Ĳ�ֵ�������˷ֱ���
% C_R_forewing=2.58; C_R_hindwing=0.9; 
% C_T=1.2;                 % 2005-JFM-Andersen_b
% A=1.4;   B=1.0;      % 2005-JFM-Andersen_b
% C_L_Andersen=C_T*sin(2*alpha_4*pi/180); 
% C_D_Andersen=A-B*cos(2*alpha_4*pi/180);
A=1.2; B=1.4;   C=1.0;
C_L_wang=A*sin(2*alpha_4*pi/180);       % �ɶ�άCFD����������϶���,
C_D_wang=B-C*cos(2*alpha_4*pi/180);  % �ɶ�άCFD����������϶���,�������������������̬����Ŀ�Ƭ
% figure(6)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_wang,'r-',alpha_4,C_D_wang,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) 2005-JFM-Andersen_a
% C_R=1.0;               % C_R��(1.0,1.4)  % 2005-JFM-Andersen_a
% C_T=1.0;                   % 2005-JFM-Andersen_a
% C_d_alpha0=0.13;              % 2005-JFM-Andersen_a  % C_d0��(0.1,0.3)
% C_d_alphapi_2=2.0;               % 2005-JFM-Andersen_a
C_T=1.833;                  % 2007-JFM-Berman  % Wang ZJ���1999-Science-Dickinson����
C_d_alpha0=0.21;                  % 2007-JFM-Berman % Wang ZJ���1999-Science-Dickinson����
C_d_alphapi_2=3.35;                 % 2007-JFM-Berman % Wang ZJ���1999-Science-Dickinson����  
C_L_andersen=C_T*sin(2*alpha_4*pi/180); 
C_D_andersen=C_d_alpha0*(cos(alpha_4*pi/180)).^2+C_d_alphapi_2*(sin(alpha_4*pi/180)).^2;  % ��������崵�����ֽ��
% C_L_andersen=C_D_andersen./tan(alpha_4*pi/180);  % ����ϵ��;�ο�2014-JRSI-Nabawy
% figure(7)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_andersen,'r-',alpha_4,C_D_andersen,'b-') 
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) 2010-JFM-Whitney
% Here mentioned coefficients is for the situation of the Reynolds number (��100), 
% geometry and flapping kinematics of a typical Drosophila wing, Dickinson et al. (1999)
% C_R=1.55;   % ת������������ϵ������ת�����ʺ�ת�����λ�ñ仯 
C_Lmax=1.8;   C_Dmax=3.4;     C_D0=0.4;  
C_L_whitney=C_Lmax*sin(2*alpha_4*pi/180);  
C_D_whitney=(C_Dmax+C_D0)/2-(C_Dmax-C_D0)/2*cos(2*alpha_4*pi/180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) 2006-IEEE TR-Deng Xinyan & 2010-EAE-RJ Wood
% C_rot_theo=2*pi*(0.75-x_0nd(r));��  % ת������������ϵ��
% C_N_deng=3.4*sin(abs(alpha_4*pi/180));     % 2010-EAE-RJ Wood           % ����1999-Science-Dickinson��ʵ������
% C_T_deng=0.4*(cos(2*alpha_4*pi/180)).^2.*(alpha_4>-45 & alpha_4<45);   % ��alpha��(-pi/4,pi/4), or else C_T=0;
C_N_deng=3.4*sin(alpha_4*pi/180);         % 2006-IEEE TR-Deng Xinyan % ����1999-Science-Dickinson��ʵ������
C_T_deng=0.4*(cos(2*alpha_4*pi/180)).^2.*(alpha_4>=0 & alpha_4<=45);   % ��alpha��[0,pi/4], or else C_T=0;
C_L_deng=C_N_deng.*cos(alpha_4*pi/180)-C_T_deng.*sin(alpha_4*pi/180);
C_D_deng=C_N_deng.*sin(alpha_4*pi/180)+C_T_deng.*cos(alpha_4*pi/180);
% figure(9)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_N_deng,'k-',alpha_4,C_T_deng,'r-',alpha_4,C_L_deng,'g-',alpha_4,C_D_deng,'b-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) 2010-DSCC-Katie_Byl
% �򻯵�Ŀ���ǣ�
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
%% (8) ������ɭ����ǰԵ�������ⷨ����������ϵ������ģ��
% C_L_polhamus=K_p*sin(alpha)*(cos(alpha)).^2+K_v*(cos(alpha))*sin(alpha).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) 2014-JRSI-Nabawy
C_La2d=2*pi;    % ��λ��rad^(-1), ����2*pi (rad^(-1))= 0.11 (deg^(-1));   
% C_La2d=5.15662;  % ��λ��rad^(-1), ��Ե���������ŵ����ƽ��, ����ʹ�� 0.09 (deg^(-1))
% The parameter E is the edge correction proposed by Jones [42] for the lifting line theory 
% and is evaluated as the quotient of the wing semi-perimeter to its length[42,43].
% lambda=0.627; % ��Ϊ����һ���������ӣ��Ա���ʵ����Զ��õ�������ϵ�����жԱ�
lambda=0.755;
E=lambda*1.146071666;   % E=C/2/R;
% The parameter k is the so-called 'k-factor' included to correct for the difference in efficiency between the assumed ideal
% uniform downwash distribution and real downwash distribution [35,40,44]. For this work, the k-factor required within
% equation (2.11) will be estimated using the induced power factor expression of hovering actuator disc models [45].
k=1.51;  % ��Թ�Ӭ
% for the span of a single  wing 
AR=3.40158;     % ע�������չ�ұ�AR��(3,5)�ǱȽϺõġ�% ��Թ�Ӭ
C_La_Ny=C_La2d/(E+k*C_La2d/(pi*AR));     % �����������������ۻ�õ���ά�����������б��  % ���: C_La_Ny =2.7506;
% C_L_nabawy=(0.5*C_La2d/(E+k*C_La2d/(pi*AR)))*sin(2*alpha_4);
C_L_nabawy=0.5*C_La_Ny*sin(2*alpha_4*pi/180);   % ��2009-JEB-Rotational accelerations ��Fig.7ͼʵ�����ݶԱȷ��Ϻܺ�
% C_D_nabawy=C_Lnabawy.*tan(alpha_4);   % C_Dnabawy=2*C_T*(sin(alpha))^2=C_La_Ny*(sin(alpha))^2;  
C_D_nabawy=C_La_Ny*(sin(alpha_4*pi/180)).^2;  
% �����ǲ���ǰԵ�������ⷨ�����Ĳ�����ɭ����ϵ��ģ��
% C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha)).*(cos(alpha)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha));
C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha_4*pi/180)); 
% figure(11)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_nabawy,'r-',alpha_4,C_D_nabawy,'b-',alpha_4,C_L_polhamus1,'g-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) 2014-AST-Taha
% a0 is the lift curve slope of the two-dimensional airfoil section, e.g. it is equal to 2�� for a flat plate or a very thin cambered shape.
% For conventional airfoils,it could be determined from lift curves such as the ones presented by Abbott and Doenhoff[1].
a_0=2*pi;
%  AR being based on one wing: AR=R^2/S;
AR=2.8;  % ע�������չ�ұ�AR��3�ǱȽϺõġ� 
C_La_Ta=pi*AR/(1+sqrt((pi*AR/a_0)^2+1));     % ����չ�����߻�õ���ά�����������б��   % ���:C_La_Ta =3.2334;
% C_L_taha=(0.5*pi*AR/(1+sqrt((pi*AR/a_0)^2+1))).*sin(2*alpha_4*pi/180);
C_L_taha=0.5*C_La_Ta*sin(2*alpha_4*pi/180);
C_D_taha=C_L_taha.*tan(alpha_4*pi/180);
% �����ǲ���ǰԵ�������ⷨ�����Ĳ�����ɭ����ϵ��ģ��
% % C_L_polhamus2=C_La_Ta*sin(alpha).*(cos(alpha)).^2+(C_La_Ta-C_La_Ta^2/(pi*AR)).*cos(alpha)*(sin(alpha)).^2;  
% C_L_polhamus2=C_La_Ta*sin(alpha_4*pi/180).*(cos(alpha_4*pi/180)).^2+(C_La_Ta-C_La_Ta^2/(pi*AR))*cos(alpha_4*pi/180).*(sin(alpha_4*pi/180)).^2; 
C_L_polhamus2=(0.5*C_La_Ta*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-C_La_Ta/(pi*AR)).*sin(alpha_4*pi/180));
cftool  % sum of sine г��������� or Fourier ������� ƽ��������ϵ��
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
%% ��һ���֡�������ϵ���ĶԱȷ���
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
%% �ڶ����֡�������ϵ���ĶԱȷ���
figure(14)
C_D=plot(alpha_4,C_D_dickinson,'k-',alpha_4,C_D_khan,'r--',...
                alpha_4,C_D_wang,'r-.',alpha_4,C_D_andersen,'r:',alpha_4,C_D_whitney,'b--',...
                alpha_4,C_D_deng,'b-.',alpha_4,C_D_byl,'b:',alpha_4,C_D_nabawy,'g--',alpha_4,C_D_taha,'g-.');
xlabel('\it\alpha (deg.)')
ylabel('\it C_D (\alpha )')
title('Aerodynamic coefficients of drag \itvs. \alpha \rm for flapping wing')
legend('\itC_{D,dickinson}','\itC_{D,khan}','\itC_{D,wang}','\itC_{D,andersen}','\itC_{D,whitney}',...
              '\itC_{D,deng}','\itC_{D,byl}','\itC_{D,nabawy}','\itC_{D,taha}')
set(C_D,'LineWidth',2)     %�μ�������MATLABԭ���빤��Ӧ�á�(�ڶ���)�߻�����/�� P182-183
grid on
%% ��������ʾ������趨
v_axis=axis;     %axis([xmin,xmax,ymin,ymax])
%Result: v_axis =  0   100     0     4
% v_axis(1)=0;          %ָ��x�����Сֵ
% v_axis(2)=100;      %ָ��x������ֵ
v_axis(3)=-0.5;          %ָ��y�����Сֵ
v_axis(4)=3.5;           %ָ��y������ֵ
axis(v_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���alpha�����ֵ����Сֵ����λ������
% [alpha_max, locat_max]=max(alpha);
% [alpha_min, locat_min]=min(alpha);
% %Result:
% % alpha_max =134.9996
% % locat_max = 11
% % alpha_min =45.0088
% % locat_min =51
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % alpha=(beta*sin(2*pi*f*t+delta))*180/pi                      %���ι��ǵ�alpha�仯����
%%% Result: alpha = 45*sin(pi/4 + (pi*t)/2)
%%%%%%%%%%%%%%%%%%%%
%���alpha�����ֵ����Сֵ����λ������
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
%���alpha�����ֵ����Сֵ����λ������
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


