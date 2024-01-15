% 2D impact oscillator with friction [main testing file]
clear;

F_scheme = 'EXP';   % explicit (EXP) or implicit (IMP) scheme
F_exc = 1;          % switch between two excitation signal sets (1 or 2)
F_kap = 1;          % cap kappa (1) to the maximum allowed value to avoid poor accuracy or not (0)
F_mc = 1;           % apply modal correction (1) or not (0)

%% TIME CONSTANTS %%%%%%%%%%%%%
dur = 0.02;

%% PHYSICAL CONSTANTS %%%%%%%%%%%%%
m = 0.001;
k = 1000;
r = 0.002;
xb = 0.0000;
if F_exc  == 1
    xb = 0.00002;
end
kap = 1e7;
alp = 1.25;
cr = 0.1;
vim = 0.2;
chi = 1e5;
thetd = 1.0;

%% PACK THE PARAMETERS INTO A SINGLE STRUCT %%%
par.m = m;
par.k = k;
par.r = r;
par.xb = xb;
par.kap = kap;
par.alp = alp;
par.cr = cr;
par.chi = chi;
par.vim = vim;
par.thetd = thetd;

%% INPUT PARAMETERS %%%%%
if F_exc == 1
    om0 = sqrt(k/m);
    fd = 1.0*om0/(2*pi);      
    par.fdx = fd;
    par.fdy = fd;
    par.ampx = -0.1;
    par.ampy = -0.1;
    par.phax = 0;
    par.phay =  0.5*pi;
elseif F_exc == 2
    par.fdx = 0;
    par.fdy = 100;
    par.ampx = 1;
    par.ampy = -1.1;
    par.phax = 0;
    par.phay =  0;
end

%% SIMULATION(S) %%%%%%%%%%%%%%%%%%%%%
% compute 'exact solution' with high sampling frequency
FsES = 256*44100;
dtES = 1/FsES;
NsES = ceil(dur*FsES);
tES = (0:(NsES-1))*dtES;
inpES = geninp(par,tES);
if F_scheme == 'IMP'
    outpES = simIMPfunc(inpES,par,FsES,F_kap,F_mc);   
else
    outpES = simEXPfunc(inpES,par,FsES,F_kap,F_mc);
end

% compute solution at audio rate
Fs = 1*44100;
dt = 1/Fs;
Ns = ceil(dur*Fs);
t = (0:(Ns-1))*dt;
inp = geninp(par,t);
if F_scheme == 'IMP'
    outp = simIMPfunc(inp,par,Fs,F_kap,F_mc);
else
    outp = simEXPfunc(inp,par,Fs,F_kap,F_mc);
end

%% PLOTTING %%%%%%%
LW = 1.5;
HF = figure(1);
clf;
HF.Position = [100,200,1800,700];

HS0 = subplot(1,2,1);
aa = 1;
HP = patch([1000*xb aa aa 1000*xb],0.995*[-aa -aa aa aa],0.9*[1 1 1]);
set(HP,'EdgeColor','None');
hold on;
plot(1000*outp.x,1000*outp.y,'.');
hold off;
grid;
axis equal;
axis([-aa aa -aa aa]);
set(gca,'Box','on');
set(gca,'Layer','top')
xlabel('x (mm)');
ylabel('y (mm)');

HS1 = subplot(4,2,2);
plot(1000*tES,1000*outpES.x,'k-','LineWidth',LW);
hold on;
plot(1000*tES,1000*outpES.y,'k-','LineWidth',LW);
h3 = plot(1000*t,1000*outp.x,'b--','LineWidth',LW);
h4 = plot(1000*t,1000*outp.y,'r--','LineWidth',LW);
hold off;
grid;
ylabel('displacement (mm)');
legend([h3 h4],'x','y','Location','NorthWest');

HS2 = subplot(4,2,4);
plot(1000*tES,outpES.Fc,'k-','LineWidth',LW);
hold on;
plot(1000*t,outp.Fc,'--','color',[0 0.8 0],'LineWidth',LW);
hold off;
grid;
ylabel('F_{c} (N)');
%axis([0 1000*dur -10 2]);

HS3 = subplot(4,2,6);
plot(1000*tES,outpES.Qc,'k-','LineWidth',LW);
hold on;
plot(1000*t,outp.Qc,'--','color',[0 0.8 0],'LineWidth',LW);
hold off;
xlabel('time (ms)');
ylabel('Q_{c} (W)');
grid;

HS4 = subplot(4,2,8);
plot(1000*tES,outpES.H,'k-','LineWidth',LW);
hold on;
plot(1000*t,outp.H,'--','color',[0 0.8 0],'LineWidth',LW);
hold off;
xlabel('time (ms)');
ylabel('H (J)');
grid;

linkaxes([HS1 HS2 HS3 HS4],'x');