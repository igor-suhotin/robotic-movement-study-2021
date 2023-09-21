global myu k1 k2 Xtol Vtol;

myu = 2.5;
k2 = 0.8;
b = 0.75;
k1 = b*k2;

filename = 'AnimationFor_k2_08.gif';

%     clr = [1 0.55 0.8]; %розовый
 %    clr = [1 1 0]; %желтый
%     clr = [0.69 0.87 0.9]; %голубой
         clr = [0 0.5 1]; %синий
% 
 %    clr0=[1 0 0.5]; %розовый
    clr0=[0 1 0.5]; %зеленый
%      clr0=[0 0.5 1]; %синий
%    clr0=[1 1 1]; %белый

Vtol=0.0000000000001;
Xtol=0.0000000001;
y0 = [0,0]; % Ќачальные услови€ (задача  оши) - вектор-столбец из величин s, fi, ds/dt, dfi/dt.
t0 = 0; % Ќачальное врем€ интегрировани€ (сек)
% Ўаг выдачи результатов интегрировани€ (сек)
n=100;
tfin=12*pi;
grstep = tfin/n;
tout=t0:grstep:tfin;
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t,y] = ode45(@MY_equation,tout,y0,opts);

% figure
% subplot(1,2,1)
% plot(t,y(:,1));
% title('X');
% subplot(1,2,2)
% plot(t,y(:,2));
% title('V');

x=y(:,1);
r=2; omega=1; R=2/3*r;
a=3*r; c=2*r;
fi=omega*t;
xO=x; yO=c/2;
xB=-a/2+x; yB=c;    xC=a/2+x; yC=yB;
xA=xB; yA=0;        xD=xC; yD=yA;
xM=xO+R*sin(fi); yM=yO-R*cos(fi);


h=figure('position',[0 0 800 300]);
axis equal;
xlim('manual');
ylim('manual');
xlim([min(x)-a max(x)+a]); % задание размера оси x
ylim([-1 c+1]); % задание размера оси y

hold on
k=linspace(0,2*pi,300); %массив размера 300 с начальным и конечным значени€ми 0 и 2*pi
circO=plot(R*cos(k)+xO(1),R*sin(k)+yO); %окружность
AB=line([xA(1) xB(1)],[yA yB], 'Color', 'Black');
BC=line([xB(1) xC(1)],[yB yC], 'Color', 'Black');
CD=line([xC(1) xD(1)],[yC yD], 'Color', 'Black');
% DA=line([xD(1) xA(1)],[yD+0.1 yA+0.1], 'Color', [0.69 0.87 0.9],'linewidth',5);
DA=line([xD(1) xA(1)],[yD+0.1 yA+0.1], 'Color', [0 1 0],'linewidth',5);
O=plot(xO(1), yO, 'ro', 'MarkerSize', 1, 'MarkerFaceColor', 'black');
% M=plot(xM(1), yM(1), 'ro', 'MarkerSize', 15,'MarkerFaceColor',clr0,'color',clr,'linewidth',5);
M=plot(xM(1), yM(1), 'ro', 'MarkerSize', 12,'MarkerFaceColor',clr0,'color',clr,'linewidth',3);
xx=line([min(x)-a max(x)+a],[0 0], 'Color', 'Black');
xy=line([-a/2 -a/2],[-1 c+1], 'Color', 'Black');
str=line([x(1)-c/2 x(1)+c/2],[-0.4 -0.4], 'Color', 'Black');
strv=line([x(1)-c/2 x(1)+c/2],[-0.4 -0.4], 'Color', 'Black');
strn=line([x(1)-c/2 x(1)+c/2],[-0.4 -0.4], 'Color', 'Black');

for i=1:length(t)
    
    set(circO,'Xdata',R*cos(k)+xO(i));
    set(AB, 'XData', [xA(i) xB(i)]);
    set(BC, 'XData', [xB(i) xC(i)]);
    set(CD, 'XData', [xC(i) xD(i)]);
    
    if abs(y(i,2)) < Vtol
        clrDA = [1 0.55 0.8];%розовый при покое
        xb=x(i)-c/2; yb=-0.4;
        xev=x(i); yev=yb;
        xen=x(i); yen=yb;
    elseif y(i,2) > 0
        clrDA = [1, 165.0/255, 0]; %оранжевый при движении вправо
        xb=x(i)+c/2; yb=-0.4;
        xev=x(i)-c/16; yev=yb+0.3;
        xen=x(i)-c/16; yen=yb-0.3;
    elseif y(i,2) < 0
        clrDA = [0.69 0.87 0.9];%голубой при движении влево
        xb=x(i)-c/2; yb=-0.4;
        xev=x(i)+c/16; yev=yb+0.3;
        xen=x(i)+c/16; yen=yb-0.3;
    end
    set(DA, 'XData', [xD(i) xA(i)],'Color',clrDA);
    set(O, 'XData', xO(i));
    set(M, 'XData', xM(i), 'YData', yM(i),'MarkerFaceColor',clr0,'color',clr);
    
    set(str,'XData', [x(i)-c/2 x(i)+c/2]);
    set(strv,'XData', [xb xev],'YData',[yb yev]);
    set(strn,'XData', [xb xen],'YData',[yb yen]);
   % pause(0.1);
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',Vtol);
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.000005);
    else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',Vtol);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.000005);
    end
end
