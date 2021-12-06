function [out]=WILCOUmodelRangeF(x,JZ,JZC,omega,LF1,UF1,NN)
% Whittle Likelihood for elliptical OU over a range of frequencies -
% fit to power spectrum and complementary spectrum
SZZ = ((1/x(3)+x(3))^2/4)*x(4)./(x(1)^2+(omega-x(2)).^2) + ((1/x(3)-x(3))^2/4)*x(4)./(x(1)^2+(omega+x(2)).^2);
RZZ = (x(4)/4)*(1/x(3)^2-x(3)^2)*(1./(x(1)^2+(omega-x(2)).^2) + 1./(x(1)^2+(omega+x(2)).^2))*exp(1i*2*x(5));

LF2 = NN-UF1+1;
UF2 = NN-LF1+1;

dets = SZZ(LF1:UF1).*SZZ(UF2:-1:LF2)-RZZ(LF1:UF1).*conj(RZZ(LF1:UF1));

out=0;
for ii = LF1:UF1
out = out+log(dets(ii+1-LF1))+[conj(JZ(ii)) conj(JZC(ii))]*[SZZ(UF2+LF1-ii) -RZZ(ii); -conj(RZZ(ii)) SZZ(ii)]*[JZ(ii); JZC(ii)]/dets(ii+1-LF1);
end