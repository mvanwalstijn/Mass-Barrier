function[inp] = geninp(par,t,F_hann)

if nargin == 2
    F_hann = 0;
end

inp.x = par.ampx*cos(2*pi*par.fdx*t + par.phax);
inp.y = par.ampy*cos(2*pi*par.fdy*t + par.phay);

if F_hann
    N = length(t);
    win = (hann(N))';
    inp.x = win.*inp.x;
    inp.y = win.*inp.y;
end




%% INPUT %%%%%%%%%%%%%%%%
% inp.x = zeros(1,Ns);
% inp.y = zeros(1,Ns);
% if F_exc == 1
%     om0 = sqrt(k/m);
%     fd = 1.0*om0/(2*pi);       % drive at resonance
%     amp = 0.1;        % amplitude
%     PO = 0.5*pi;      % Phase Offset
%     inp.x = -amp*cos(2*pi*fd*t);
%     inp.y = -amp*cos(2*pi*fd*t + PO);
%     %inp.y = zeros(size(inp.x));
% else
%     fd = 100;       % drive at resonance
%     amp = 1.1;        % amplitude
%     inp.x = ones(1,Ns);
%     % Nf = ceil(0.01*dur*Fs);
%     % win =  [(0:(Nf-1))/Nf ones(1,Ns-Nf)];
%     inp.y = -amp*sin(2*pi*fd*t);
% end

