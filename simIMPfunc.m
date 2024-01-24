function[outp] = simIMPfunc(inp,par,Fs,F_kap,F_mc)
% 2D impact oscillator with friction [implicit scheme, function]

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
    fprintf(1,'\n IMP: kappa = %1.1f capped to %1.1f',kap, kap_max);
    kap = kap_max;
end

%% INITIAL VALUES %%%%%%%%%
u = [0; 0];
um = [0; 0];

%% INPUT & OUTPUT %%%%%%%%%%%%%%%%
outp.t = t;
outp.x = zeros(1,Ns);
outp.y = zeros(1,Ns);
outp.Fc = zeros(1,Ns);
outp.Fr = zeros(1,Ns);
outp.h = zeros(1,Ns);
outp.v = zeros(1,Ns);
outp.iter = zeros(1,Ns);
outp.H = zeros(1,Ns);
outp.Ff = zeros(1,Ns);

%% SOLUTION CONSTANTS %%%%%%%%%%%%
EPS = eps;
EPS2 = eps^2;
maxiter = 1000;
TH = 100*eps;
s = 0.000001*[1; 1];

%% TIME-STEPPING LOOP %%%%%%%%%%%%%
eta = 1;
hm = eta;

for n=1:Ns
    %required variables 
    fe = [inp.x(n); inp.y(n)];
    z = 0.5*(B*um - A*u - C*fe);
    us = um - 2*z;
    vi = (us(1) - u(1))/dt;         % estimate velocity between n and n+1
    if us(1) > xb && u(1) <= xb      % when making contact, set eta
        eta = ((0.2*alp + 1.3)*(1-cr))/((cr + chi*dt^2)*(vi + chi*dt^2));  % adjusted Sun et al expression
    end
    %x component
    wstar = 0.5*(u(1) + um(1)) - xb;               % past contact variable    
    Phimh = (kap/(alp+1))*max(0,wstar)^(alp+1);
    dPhimh = kap*max(0,wstar)^alp;
    ddPhimh = kap*alp*max(0,wstar)^(alp-1);
    zx = z(1);
    sx = s(1);
    diff = 1;
    iter = 0;
    while diff > TH & iter < maxiter
        Phiph = (kap/(alp+1))*max(0,0.5*sx+wstar)^(alp+1);
        dPhiph = kap*max(0,0.5*sx+wstar)^alp;                   
        v = sx/(2*dt);
        if v < -1/eta
            h = -1/v;   
        else
            h = eta;
        end
        % dh = (0.5*bet^2*dt)/(dt^2 + (0.5*bet*sx)^2);      % analytic dh
        dh = (h-hm)/dt;     % Calculate dh numerically instead 
        hm = h;
        q =  C*(1 + (0.5*sx/dt)*h);
        dq = (0.5*C/dt)*(h + sx*dh);
        Fcx = 2*(((Phiph - Phimh)*sx)/(sx^2 + EPS2) + 0.5*dPhimh*(EPS2/(sx^2 + EPS2)));
        dFcx = 2*((0.5*dPhiph*sx - Phiph + Phimh)/(sx^2 + EPS2) + (EPS2/(sx^2 + EPS2))*(0.125*ddPhimh));
        Fnl = sx + 2*zx + q*Fcx;
        dFnl = 1 + q*dFcx + dq*Fcx;
        sxnew = sx - Fnl/dFnl;
        diff = abs(sxnew - sx);
        sx = sxnew;
        iter = iter + 1;
    end
    outp.iter(n) = iter;

    if iter >= maxiter
        disp('maxiter reached');
    end

    Phiph = (kap/(alp+1))*max(0,0.5*sx+wstar)^(alp+1);  
    % Calculate Fr using a hard switch for abs(sx) <= eps
    if abs(sx) < eps
        Fr = -dPhimh;
    else
        Fr = -2*(Phiph - Phimh)/sx;
    end  
    rhox = (sx/(2*dt))*h;
    s(1) = sx;
    
    %update the y component
    v = sx/(2*dt);
    if v < -1/eta
        h = -1/v;   
    else
        h = eta;
    end
    rhoy = (sx/(2*dt))*h;
    Fcy = Fr*(1 + rhoy);
    zy = z(2);
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
    s(2) = sy;
    
    %update state variables 
    up = s + um;

    %record some outputs
    outp.x(n) = u(1);
    outp.y(n) = u(2);
    outp.Fr(n) = Fr;
    outp.Fc(n) = (1 + rhox)*Fr;
    outp.rhox(n) = rhox;
    outp.Ff(n) = thet*(1 + rhoy)*Fr;
    outp.q(n) = q;
    outp.h(n) = h;
    outp.vx(n) = sx/(2*dt);
    outp.vy(n) = sy/(2*dt);
    outp.Qc(n) = -( ((sx/(2*dt))^2)*h + (sy/(2*dt))*(1 + rhoy)*thet )*Fr;
    outp.H(n) = 0.5*m*(up - u)'*(up - u)/dt^2 + 0.5*khat*(up + u)'*(up + u)/4 + Phiph;
    
    %memorise
    um = u;
    u = up;
end


