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
outp.qy = zeros(1,Ns);
outp.h = zeros(1,Ns);
outp.vx = zeros(1,Ns);
outp.vy = zeros(1,Ns);
outp.Qc = zeros(1,Ns);
outp.H = zeros(1,Ns);
outp.Ff = zeros(1,Ns);

%% SOLUTION CONSTANTS %%%%%%%%%%%%
EPSIL = eps;

%% TIME-STEPPING LOOP %%%%%%%%%%%%%
eta = 1;
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
        h = -1/vh;
    else
        h = eta;
    end
    gam = (0.5/dt)*h;
    lam = zx - 0.5*gam*C*psimh^2;    
    
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
    a = 0.25*C*gam*g^2;
    b = 1 + C*gam*psimh*g + 0.25*C*g^2;
    c = C*psimh*g + 2*zx;
    D = b^2 - 4*a*c;
    sx = -2*c/(b + sqrt(D));
    
    %update the auxiliary variable, the restoring force, and qx;
    psiph = psimh + 0.5*g*sx;           
    Fr = -0.5*(psiph+psimh)*g; 
    rhoy = ((w - wm)/dt)*h;
    qy =  C*(1 + rhoy);
       
    %update the y component
    zy = z(2);
    p = -qy*Fr;
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

    %update state variables 
    s(1) = sx;      
    s(2) = sy;
    up = s + um;

    %record
    outp.x(n) = u(1);
    outp.y(n) = u(2);
    outp.Fc(n) = (1 + (0.5*sx/dt)*h)*Fr;
    outp.rhox(n) = gam*sx;
    outp.Ff(n) = thet*(1 + rhoy)*Fr;
    outp.Fr(n) = Fr;
    outp.qy(n) = qy;
    outp.h(n) = h;
    outp.vx(n) = sx/(2*dt);
    outp.vy(n) = sy/(2*dt);
    outp.Qc(n) = -( ((sx/(2*dt))^2)*h + (sy/(2*dt))*(1 + rhoy)*thet )*Fr;
    outp.H(n) = 0.5*m*(up - u)'*(up - u)/dt^2 + 0.5*khat*(up + u)'*(up + u)/4 + 0.5*psiph^2;

    %memorise
    um = u;
    u = up;
    psimh = psiph;
end

