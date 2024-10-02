function [Fp0, Fm0, Gp0, Gm0, Qp0, Qm0, par]=SFMSA_init(par,ind)
%% set initial value of symmetric rotating wave for SFMSA

% set up rotating wave
G0=(1+par(ind.epsa)-(par(ind.a)*par(ind.A)-par(ind.B))/(1-par(ind.a)))/2; % initial guess for G0 on ellipse
Q0=par(ind.B)*G0/((par(ind.a)*par(ind.A)+(1-par(ind.a))*G0));
F0=[sqrt((par(ind.A)-G0)./G0);0]; % in corotating frame imaginary part is zero

% set parameters such that above is rotating wave
omega0 = par(ind.alpha)*G0-par(ind.beta)*Q0-par(ind.epsp);
phi0 = omega0*par(ind.tau) + tan((par(ind.alpha)*G0-par(ind.beta)*Q0-par(ind.epsp) - omega0)./(G0-Q0-1-par(ind.epsa)));

% warning if G0 was a bad guess
if imag(par(ind.omega))~=0 || imag(par(ind.phi))~=0
    error('in SFMSA_init omega or psi not real');
end

% psi mod 2pi
par(ind.omega)= omega0; par(ind.phi)=mod(phi0,2*pi);

% set variables counterparts
Fp0=F0; Fm0=F0; Gp0=G0; Gm0=G0; Qp0=Q0; Qm0=Q0;
end