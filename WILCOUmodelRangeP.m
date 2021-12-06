function [out]=WILCOUmodelRangeP(x,SZ,omega,LF1,UF1)
% Whittle Likelihood for circular OU over a range of frequencies
SZZ = x(3)./(x(1)^2+(omega-x(2)).^2);
out = sum(log(SZZ(LF1:UF1)) + SZ(LF1:UF1)./SZZ(LF1:UF1));