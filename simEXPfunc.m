function[outp] = simEXPfunc(inp,par,Fs,F_kap,F_mc)
% 2D impact oscillator with friction [explicit scheme]

%% TIME CONSTANTS %%%%%%%%%%%%%
dt = 1/Fs;
Ns = length(inp.x);
t = (0:(Ns-1))*dt;

%% PHYSICAL CONSTANTS %%%%%%%%%%%%%
m = par.m;
k = par.k;
r = par.r;
xb = par.xb;
kap = par.kap;
alp = par.alp;
cr = par.cr;
chi = par.chi;
thetd = par.thetd;

%% SYSTEM PARAMETERS & COEFFICIENTS
om0 = sqrt(k/m);
sigm = r/(2*m);
om = sqrt(om0^2 - sigm^2);
OM = cos(om*dt);
UP = exp(-sigm*dt);
xi = (dt^2)/m;
if F_mc == 1
    A = 2*OM*UP;
    B = 1 + UP^2;
    C = 0.25*xi*(1 + UP^2 + 2*UP*OM);
    khat = (4*m/dt^2)*(B - A)/(B + A);
    rhat = (2*m/dt)*(4 - 2*B)/(B + A);
else
    a = (k*dt^2)/(4*m);
    b = (r*dt)/(2*m);
    A = (2 - 2*a)/(1 + a + b);
    B = (2 + 2*a)/(1 + a + b);
    C = xi/(1 + a + b);
    khat = k;
    rhat = r;
end

%% CAP KAPPA %%%%%%%%%%%%%%%%%
vim = par.vim;   % expected max impact velocity
zeta = (sqrt(pi)*gamma(1 + 1/(alp+1)))/(gamma(0.5 + 1/(alp+1)));
Nmin = 6;
kap_max = 0.5*(vim^(1 - alp))*(alp + 1)*m*((2*zeta)/(Nmin*dt))^(alp + 1);
if F_kap && kap > kap_max
    fprintf(1,'\n EXP: kappa = %1.1f capped to %1.1f',kap, kap_max);
    kap = kap_max;
end

%% INITIAL VALUES %%%%%%%%%
u = [0; 0];
um = [0; 0];
psimh = 0;
s = [0; 0];

%% INPUT & OUTPUT %%%%%%%%%%%%%%%%
outp.t = t;
outp.x = zeros(1,Ns);
outp.y = zeros(1,Ns);
outp.Fc = zeros(1,Ns);
outp.Fr = zeros(1,Ns);
outp.h = zeros(1,Ns);
outp.vx = zeros(1,Ns);
outp.vy = zeros(1,Ns);
outp.Qc = zeros(1,Ns);
outp.H = zeros(1,Ns);
outp.Ff = zeros(1,Ns);
outp.psi = zeros(1,Ns);

%% SOLUTION CONSTANTS %%%%%%%%%%%%
EPSIL = eps;
eta = 1;

%% TIME-STEPPING LOOP %%%%%%%%%%%%%
for n=1:Ns
    %required variables
    fe = [inp.x(n); inp.y(n)];
    z = 0.5*(B*um - A*u - C*fe);
    zx = z(1);
    w = u(1) - xb;                     % current contact variable
    wm = um(1) - xb;                   % previous contact variable    
    vh = (w - wm)/dt;
    us = um - 2*z;
    vi = (us(1) - u(1))/dt;         % estimate velocity between n and n+1
    if us(1) > xb && u(1) <= xb      % when making contact, set eta
       eta = ((0.2*alp + 1.3)*(1-cr))/((cr + chi*dt^2)*(vi + chi*dt^2));  % adjusted Sun et al expression
    end
    if vh < -1/eta
        hx = -1/vh;
    else
        hx = eta;
    end
    gam = (0.5/dt)*hx*psimh;
    lam = zx - gam*C*psimh;    
    
    %calculate allowable g \in [0, gn]
    gp = 2*(sqrt(lam^2 + C*psimh^2) - lam)/(sqrt(EPSIL*C*lam^2 + (C*psimh)^2));
    if w > 0 && kap > 0
        gn = sqrt(0.5*(alp +1)*kap*w^(alp - 1));
        g = min(gn,gp);
    elseif w <= 0 && wm > 0
        g = gp;
    else
        g = 0;
    end

    %calculate sx
    sx = -(2*zx + C*psimh*g)/(1 + C*gam*g + 0.25*C*g^2);
    
    %update the auxiliary variable, the restoring forcex;
    psiph = psimh + 0.5*g*sx;   
    %psiph = max(0,psiph);        % if needed, mitigate finite-precision round-off    
    Fr = -0.5*(psiph+psimh)*g;  % contact force    
       
    %update the y component
    zy = z(2);
    vx = sx/(2*dt);
    if vx < -1/eta
        hy = -1/vx;
    else
        hy = eta;
    end   
    Fcy = Fr*(1 + vx*hy);
    p = -C*Fcy;
    if zy < -0.5*p*thetd
        sy = -p*thetd - 2*zy;
        thet = thetd;
    elseif zy > 0.5*p*thetd
        sy = p*thetd - 2*zy;
        thet = -thetd;
    else
        sy = 0;
        thet = -2*z(2)/p;
    end
    vy = sy/(2*dt);

    %update state variables 
    s(1) = sx;      
    s(2) = sy;
    up = s + um;

    %record
    outp.x(n) = u(1);
    outp.y(n) = u(2);
    outp.Fc(n) = Fcy;
    outp.Ff(n) = thet*Fcy;
    outp.Fr(n) = Fr;   
    outp.h(n) = hx;
    outp.psi(n) = psiph;
    outp.vx(n) = vx;
    outp.vy(n) = vy;
    outp.Qc(n) = vx*((psiph-psimh)/dt)*hx*psimh + ((psiph+psimh)/2)*g*(1 + vx*hy)*vy*thet;
    outp.H(n) = 0.5*m*(up - u)'*(up - u)/dt^2 + 0.5*khat*(up + u)'*(up + u)/4 + 0.5*psiph^2;

    %memorise
    um = u;
    u = up;
    psimh = psiph;
end

