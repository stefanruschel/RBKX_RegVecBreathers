function y = SFMSA(Fp, Fptau, Fm, Fmtau, Gp, Gm, Qp, Qm, p, ind)
%% r.h.s of modified spin flip model with saturable absorber and delay
%
Fpdot = 1./2*(( (1+1i*p(ind.alpha)).*Gp - (1+1i*p(ind.beta)).*Qp - 1).*Fp ...
         - (p(ind.epsa)+1i*p(ind.epsp)).*Fm + p(ind.cp)*exp(1i*p(ind.phi)).*Fptau);
Fmdot = 1./2*(( (1+1i*p(ind.alpha)).*Gm - (1+1i*p(ind.beta)).*Qm - 1).*Fm ...
         - (p(ind.epsa)+1i*p(ind.epsp)).*Fp + p(ind.cp)*exp(1i*p(ind.phi)).*Fmtau);
Gpdot = p(ind.gammaG)*( p(ind.A) - Gp.*(1+abs(Fp).^2) ) - p(ind.gammaGSF).*(Gp-Gm);
Gmdot = p(ind.gammaG)*( p(ind.A) - Gm.*(1+abs(Fm).^2) ) - p(ind.gammaGSF).*(Gm-Gp);
Qpdot = p(ind.gammaQ)*( p(ind.B) - Qp.*(1+p(ind.a).*abs(Fp).^2) ) - p(ind.gammaQSF).*(Qp-Qm);
Qmdot = p(ind.gammaQ)*( p(ind.B) - Qm.*(1+p(ind.a).*abs(Fm).^2) ) - p(ind.gammaQSF).*(Qm-Qp);     
% output
y=cat(1,real(Fpdot),imag(Fpdot),real(Fmdot),imag(Fmdot),real(Gpdot),real(Gmdot),real(Qpdot),real(Qmdot));
end
