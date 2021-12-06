function [out]=WILCOUmodelFull2(x,SZ,omega)
% Whittle Likelihood for elliptical OU over all frequencies -
% fit to power spectrum only, including aliasing effects
K=10;
omegaM = zeros(2*K+1,length(omega));
for kk = 1:K
omegaM(K+1-kk,:) = omega-2*pi*kk;
omegaM(K+1+kk,:) = omega+2*pi*kk;
end
omegaM(K+1,:) = omega;


SZZ = sum(((1/x(3)+x(3))^2/4)*x(4)./(x(1)^2+(omegaM-x(2)).^2) + ((1/x(3)-x(3))^2/4)*x(4)./(x(1)^2+(omegaM+x(2)).^2));


out = sum(log(SZZ) + SZ./SZZ);