function [out]=WILCOUmodelFullF(x,JZ,JZC,omega,NN)
% Whittle Likelihood for elliptical OU over all frequencies -
% fit to power spectrum and complementary spectrum, including aliasing
% effects
K=10;
omegaM = zeros(2*K+1,length(omega));
for kk = 1:K
omegaM(K+1-kk,:) = omega-2*pi*kk;
omegaM(K+1+kk,:) = omega+2*pi*kk;
end
omegaM(K+1,:) = omega;


SZZ = sum(((1/x(3)+x(3))^2/4)*x(4)./(x(1)^2+(omegaM-x(2)).^2) + ((1/x(3)-x(3))^2/4)*x(4)./(x(1)^2+(omegaM+x(2)).^2));
RZZ = sum((x(4)/4)*(1/x(3)^2-x(3)^2)*(1./(x(1)^2+(omegaM-x(2)).^2) + 1./(x(1)^2+(omegaM+x(2)).^2))*exp(1i*2*x(5)));

LIM = floor(NN/2);

dets = SZZ(1:(LIM+1)).*SZZ(NN:-1:(LIM+1))-RZZ(1:(LIM+1)).*conj(RZZ(1:(LIM+1)));

out=0;
for ii = 1:LIM
out = out+log(dets(ii))+[conj(JZ(ii)) conj(JZC(ii))]*[SZZ(NN+1-ii) -RZZ(ii); -conj(RZZ(ii)) SZZ(ii)]*[JZ(ii); JZC(ii)]/dets(ii);
end
out = out+0.5*(log(dets(LIM+1))+[conj(JZ(LIM+1)) conj(JZC(LIM+1))]*[SZZ(LIM+1) -RZZ(LIM+1); -conj(RZZ(LIM+1)) SZZ(LIM+1)]*[JZ(LIM+1); JZC(LIM+1)]/dets(LIM+1));
out = abs(out);