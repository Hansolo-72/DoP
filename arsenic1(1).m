function f=arsenic1 (ini)

load('expt.mat','exp')
experi=exp(:,end);                                % y variable
c0_inlet=5000;                                  % (initial conc in ppm)
L=8;                                       %input('enter the length of column (m) : ');
Q_inlet=4.2;                                %input('enter the flowrate (LPD) : ');

Q=Q_inlet*(10^-3)/(24*3600);    % flow rate in m3/sec
c0=c0_inlet*(10^-3);            % inlet conc in Kg/m3

tim = 3000;                       % time in hours

% do not change Dl
%Dl=3*10^-11;%input('enter the axial dispersion coeff. (m2/s) : ');

%Dp=1*10^-10;%input('enter the pore diffusion coeff. m2/s) : ');
%Kf=1.2*10^-6;%input(enter the mass transfer coeff. (m/s):');


Rb=0.20/2;                          % radius of bed (m)
ep=0.4;                            % particle (adsorbent) porosity
eb=0.6;                            % bed porosity
Rp=350*(10^-6);                     % particle radius (m)
v=Q/(3.1416*Rb*Rb*eb);              % velocity (m/s)

tou = tim*v*3600/L;                 % normalised time
t_n=1500;  % no of divisions of total time
l_n=10;
cpr=ones(t_n,l_n);
cpf=ones(t_n,l_n);




%-------------------------------------------------------------------------

%      Kf=ini(1);%*(10^-6);
%       Dp=ini(2);%*(10^-10);
%       Dl=ini(3);%*(10^-11);

       Kf=1.5*10^-6;
       Dp=8*10^-10;
       Dl=4*(10^-11);


Pe=v*L/Dl;                          % Peclet number
Bi=Kf*Rp/(ep*Dp);                   % Biot number
neta=ep*Dp*L/(Rp*Rp*v);             % dimensionless group neta
zi=3*Bi*neta*(1-eb)/eb;             % dimensionless group zi


    cpr=cpf;
    m=1;
    x=linspace(0,1,l_n);
    t=linspace(0,tou,t_n);
    u3=pdepe(m,@system2,@initial1,@bc1,x,t);
    cb=u3(:,end,end);
    m1=2;
    for i=1:1:l_n
    cp=pdepe(m1,@system3,@initial2,@bc2,x,t);
    cpf_1=cp(:,:,end);
    cpf(:,i)=cpf_1(:,end);
    end
      f=0;
      for y=1:1:37
          f=f+(((experi(y)-cpf(y))/cpf(y))^2);
     end
%     display(error1);


figure
%plot(t_n,f,'O')
hold on
plot(t*L/(v*3600),cb(:,end))
title('plot of conc');
xlabel('time (hours)');
ylabel('C/C0');
ylim([0 1]);
xlim([0 100]);

%--------------------------------------------------------------------------


function [c,b,s]=system2(x,t,u,DuDx)
c=1;
b=DuDx/Pe;
s=-(DuDx*(1+(1/(x*Pe)))+zi*(u-cpr(round(((t/tou))*(t_n-1)+0.5),round(((x/1))*(l_n-1)+0.5))));
end

  
function [p1,q1,pr,qr]=bc1(x1,u1,xr,ur,t)
p1=Pe*(1-u1);
q1=Pe;
pr=0;
qr=Pe;
end


function value=initial1(x)
if x < 1e-04
value=1;
else 
    value =0;
end
end


%--------------------------------------------------------------------------


function [c,b,s]=system3(x,t,u,DuDx)
d=1050;
ep=0.37;
qm=21.6*(10^-3);
k=20.6*(10^3);
%c0=arrayfun(feed,data);
c=(ep+(1-ep)*((d*qm*k)/((1+k*u*c0)^2)))/neta;
b=DuDx;
s=0;
end


function [p1,q1,pr,qr]=bc2(x1,u1,xr,ur,t)
p1=0;
q1=1;
pr=-Bi*(cb(round(((t/tou))*(t_n-1)+0.5),round(((x1/1))*(l_n-1)+0.5))-ur);
qr=1;
end



function value=initial2(x)
if x < 1e-04
value=1;
else 
value =0;
end
end

       
end





    
