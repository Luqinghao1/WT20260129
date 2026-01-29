function Dp = Composite_shale_oil_reservoir_fitfun(x,t)
% Ver:2025-12-1
% -------------------
% x = [1e-3,1e-4,1000,100,0.1,4,10,0.4,0.08,1e-3];
x = [1e-3,1e-4,1000,100,0.1,4,0.4,0.08,1e-3,10];
t = logspace(-3,3,100);
clear global
% -------------------
global hDp hDp_der
% 拟合未知参数
kf = x(1);  % 1;   % mD,   1e-15;  内区的渗透率
km = x(2);  % 0.1; % mD, 0.1e-15;  外区的渗透率

L = x(3);  % 1000;
Lf = x(4);  % 100;  裂缝半厂

LfD = x(5);  % Lf/L;
rmD = x(6);  % rm/L;

omga1 = x(7);   %  0.4;
omga2 = x(8);   % 0.08;
remda1 = x(9);  % 1e-3;    %% 窜流系数，基岩到裂缝的拟稳态窜流

reD=x(10);

% 可知参数
nf = 4;   % 裂缝条数
xwD = linspace(-0.9,0.9,nf);  % 裂缝位置
phi = 0.05;   % 基质孔隙度   
h = 20;     % 有效厚度
mu = 0.5;   % 粘度
B = 1.05;   % 体积系数
Ct = 5e-4;  % 综合压缩系数
q = 5;     % 定产生产(m3/s)

M12 = kf/km;  % 流度比

%关于时间tD的Laplace逆变换   %时间t以天为单位
tD =  14.4*kf*t/(phi*mu*Ct*L^2);  % t = phi_m*mu*Ct*L^2*tD/(14.4*kf);
PD=zeros(size(tD));

N=4;
% n = length(tD);
for i=1:length(tD)
     PD(i)=0; 
     y = 0;
     
     for j=1:N
        z=(log(2)/tD(i))*j;   
        v=0;
        for k=floor((j+1)/2):min(N/2,j)
            v=v+power(k,N/2+1)*factorial(2*k)/(factorial(N/2-k)*(factorial(k))^2*factorial(j-k)*factorial(2*k-j));
        end
        if(mod(N/2+j,2)==1), v=-v; else,   v=v; end
        temp = omga2;
        fs1 = omga1+remda1*temp/(remda1+z*temp);
        fs2 = M12*temp;
        y = PWD_inf(z,fs1,fs2,M12,LfD,rmD,nf,xwD,reD);
        pf = y(nf+1);
        % ------------
        PD(i) = PD(i)+v*pf*log(2)/tD(i);
     end
     % -------------------
     % 摄动法考虑压敏
     gamaD = 0.02;
     if gamaD == 0
         PD(i) = PD(i);  else
         PD(i) = -1/gamaD*log(1-gamaD*PD(i));
     end
end

% t = phi_m*mu*Ct*L^2*tD/(14.4*kf);
Dp = 1.842e-3*q*mu*B*PD/(kf*h);    
Dp(isnan(Dp)) = 0;   
Dp_der = t(2:end).*diff(Dp)./diff(t);
if isempty(hDp)
    hDp = loglog(t,Dp,'-');  hold on
    hDp_der = loglog(t(2:end),Dp_der,'-');  
else
    set(hDp,'ydata',Dp);
    set(hDp_der,'ydata',Dp_der);
end
drawnow
% grid on
% xlabel('t,h');
% ylabel('\Deltap&d(\Deltapp)/d(ln(\Deltapp)),MPa');
% title('有因次压力响应')
% legend('压力','压力导数')

end 

% ===================================================================
function pf = PWD_inf(z,fs1,fs2,M12,LfD,rmD,nf,xwD,reD)

% xwD = linspace(-0.9,0.9,nf);
ywD = zeros(size(xwD));

gama1 = sqrt(z*fs1);
gama2 = sqrt(z*fs2);
% --------------------------------------
% mAB=0;   %% 无穷大外边界
mAB=besselk(1,gama2*reD)/besseli(1,gama2*reD);  %% 封闭边界
% mAB=-besselk(0,gama2*reD)/besseli(0,gama2*reD); %% 定压边界
% --------------------------------------
Acup = M12*gama1*besselk(1,gama1*rmD)*(mAB*besseli(0,gama2*rmD)+besselk(0,gama2*rmD))+gama2*besselk(0,gama1*rmD)*(mAB*besseli(1,gama2*rmD)-besselk(1,gama2*rmD));
Acdown = M12*gama1*besseli(1,gama1*rmD)*(mAB*besseli(0,gama2*rmD)+besselk(0,gama2*rmD))-gama2*besseli(0,gama1*rmD)*(mAB*besseli(1,gama2*rmD)-besselk(1,gama2*rmD));
Ac = Acup/Acdown;

% 压力表达式
A = zeros(nf+1,nf+1);
for i=1:nf
    for j=1:nf
        y11=@(a)besselk(0,gama1*sqrt(((xwD(i)-xwD(j)-a).^2+(ywD(i)-ywD(j)).^2)))+Ac*besseli(0,gama1*sqrt(((xwD(i)-xwD(j)-a).^2+(ywD(i)-ywD(j)).^2)));
        pfD(i,j)=quadgk(y11,-LfD,LfD)/(M12*z*2*LfD);
        A(i,j)=z*pfD(i,j);
    end
end
A(:,nf+1) = -1;
A(nf+1,:) = z;
A(nf+1,nf+1) = 0;

b = zeros(nf+1,0);
b(nf+1) = 1;
pf = (A+eps(nf))\b';

% -----------------------
% 考虑井筒储存CD和表皮S
CD = 0.01;
S = 1;
pf(nf+1)=(z*pf(nf+1)+S)/(z+CD*z^2*(z*pf(nf+1)+S));
% -------------------------
end  % function pf=PWD(z,fs1,fs2,M12)0.4