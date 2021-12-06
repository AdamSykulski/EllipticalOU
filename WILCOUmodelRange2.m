function [out]=WILCOUmodelRange2(x,SZ,omega,LF1,UF1,LF2,UF2)
% Whittle Likelihood for elliptical OU over a range of frequencies -
% semi-parametric fit to power spectrum only
SZZ = ((1/x(3)+x(3))^2/4)*x(4)./(x(1)^2+(omega-x(2)).^2) + ((1/x(3)-x(3))^2/4)*x(4)./(x(1)^2+(omega+x(2)).^2);
out = sum(log(SZZ(LF1:UF1)) + SZ(LF1:UF1)./SZZ(LF1:UF1)) + sum(log(SZZ(LF2:UF2)) + SZ(LF2:UF2)./SZZ(LF2:UF2));