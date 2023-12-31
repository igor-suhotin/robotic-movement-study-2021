global myu k1 k2 Xtol Vtol;
myu = 2.5 %1.08 %1.3; % 1.5;
k2 = 0.45 %0.5 %0.44 % 0.2; %0.6; %k1 � k2 - ������������ ������. ������������, ��� ��������� k1 � k2 ����� b.
b = 0.75; %0.75;
k1 = b*k2;
Vtol=0.0000000001;
Xtol=0.0000000001;
X0 = [0,0]; % ��������� ������� (������ ����) - ������-������� �� ������� s, fi, ds/dt, dfi/dt.
t0 = 0; % ��������� ����� �������������� (���)
n=10000;% ��� ������ ����������� �������������� (���)
tfin=4*pi;
grstep = tfin/n;
tout=t0:grstep:tfin;
opts = odeset('RelTol',1e-12,'AbsTol',1e-14); %������� ����� � ������ RelTol, ������� ����� ��������������� �������� 1e-12 
% ����� ������� ����� � ������ AbsTol ������� ����� ��������������� �������� 1e-14
tic
[t,X] = ode45(@MY_equation,tout,X0,opts); 
toc


figure

subplot(1,2,1) %���������� �������� ������������ ���� �� ��������� �������
hold on
plot(t,X(:,1));
title('X');
xlabel('�����') 
ylabel('���������� X')
%whitebg('w') 

subplot(1,2,2)
hold on
plot(t,X(:,2));
title('V');
xlabel('�����') 
ylabel('�������� �������� V')