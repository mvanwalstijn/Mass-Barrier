% 2D impact oscillator with friction [test coefficient of restitution]
clear;

F_exc = 1;      % switch between two excitation signal sets
F_kap = 1;      % cap kappa to the maximum allowed value to avoid poor accuracy
F_mc = 1;       % apply modal correction

%% TIME CONSTANTS %%%%%%%%%%%%%
dur = 0.002;

%% PHYSICAL CONSTANTS %%%%%%%%%%%%%
m = 0.001;
k = 0;
r = 0.0;
xb = 0.0001;
kap = 1e9;
alp = 1.25;
vim = 0.5;
chi = 1e7;
thetd = 0.0;

%% PACK THE PARAMETERS INTO A SINGLE STRUCT %%%
par.m = m;
par.k = k;
par.r = r;
par.xb = xb;
par.kap = kap;
par.alp = alp;
par.chi = chi;
par.vim = vim;
par.thetd = thetd;


%% SIMULATION(S) %%%%%%%%%%%%%%%%%%%%%
OF = 32;
Fs = OF*44100;
dt = 1/Fs;
Ns = ceil(dur*Fs);
t = (0:(Ns-1))*dt;
inp.x = zeros(1,Ns);
inp.y = zeros(1,Ns);
durh = 0.0002;
Nh = ceil(durh*Fs);
inp.x(1:Nh) = hanning(Nh);
%inp.x(1) = 2;
Nt = 40;
crv = (1:Nt)/Nt;
sig = zeros(Nt,Ns);
cre = zeros(1,Nt);
n1 = round(0.00075*Fs);
n2 = n1 + 1;
n3 = round(0.00150*Fs)
n4 = n3 + 1;
for i=1:Nt
    par.cr = crv(i);
    outp = simEXPfunc(inp,par,Fs,F_kap,F_mc);
    sig(i,:)  = outp.x;
    vin = (sig(i,n2) - sig(i,n1))/dt;
    vout = (sig(i,n4) - sig(i,n3))/dt;
    cre(i) = -vout/vin;
end

%% PLOTTING %%%%%%%%%%%%%
HF = figure(1);
clf;

subplot(1,2,1);
plot(1000*outp.t,1000*sig,'.-');
grid;
xlabel('time (ms)');
ylabel('x (mm)');

subplot(1,2,2);
plot(crv,crv,'k-');
hold on;
plot(crv,cre,'r.','MarkerSize',10);
hold off;
grid;
xlabel('$c_{r}$','interpreter','latex');
ylabel('$\hat{c}_{r}$','interpreter','latex');
title(['Explicit (' num2str(OF) ' x oversampling)']);