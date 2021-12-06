%% INPUTS
Delta = 1/1000; % Euler Maruyama increment
N = 1859; % Time series length (plus 100 as we will lose these points for simulation burn in)
rng(1); % random number seed
alpha1 = 0.02; beta1 = 1; alpha2 = -0.5; beta2 = -0.3; A1 = 2; % Elliptical OU parameters
TT = 1000; % number of repeats
%% PRELIMINARIES
R1 = A1*(-beta2-1i*alpha2)/beta1; % noise relation
betaP = sign(beta1)*sqrt(beta1^2-beta2^2-alpha2^2); % peak frequency
T = N/Delta; % total number of increments in Euler Maruyama generation
PP = zeros(TT,5); BTS1 = zeros(TT,5); BTS2 = zeros(TT,5); % other initial entries
for ii = 1:TT
    ii
ZZ = zeros(1,T);
for t = 1:T-1 % Euler Maruyama generation
    s1 = (A1+real(R1))/2; s2 = (A1-real(R1))/2; % magnitude of noisy in x/y components
    corr=imag(R1)/(2*sqrt(s1)*sqrt(s2)); % correlation between x/y noise components (please refer to Sykulski and Percival (2016) IEEE paper for details)
    B1 = randn(1,1); C1 = randn(1,1); D1 = corr*B1 + sqrt(1-corr^2)*C1; % correlated noise in x/y
    N1 = sqrt(s1)*B1; N2 = sqrt(s2)*D1; % scaled to correct amplitudes
    ZZ(t+1) = ZZ(t) - (alpha1-1i*beta1)*ZZ(t)*Delta - (alpha2-1i*beta2)*conj(ZZ(t))*Delta + sqrt(Delta)*(N1+1i.*N2); % Euler Maruyama forwards equation
end
ZNew3=ZZ(1/Delta:1/Delta:end); % subsample to correct sampling frequency
ZNew4 = ZNew3(101:end); NN = length(ZNew4); % throw away first 100 points for burn in
omegaN=0:2*pi/NN:2*pi*(1-1/NN); omegaN=fftshift(omegaN); omegaN(1:floor(NN/2))=omegaN(1:floor(NN/2))-2*pi;
%%
JZ = (1/sqrt(NN))*fftshift(fft(ZNew4)); % Fourier transform
JZC = (1/sqrt(NN))*fftshift(fft(conj(ZNew4))); % Fourier transform of conjugate
spec1=(1/NN)*fftshift(abs(fft(ZNew4)).^2); % periodogram
spec2=(1/NN)*fftshift((fft(ZNew4)).*conj(fft(conj(ZNew4)))); % complementary spectrum
options=optimset('gradobj','on','MaxFunEval',10000,'MaxIter',10000,'TolX',1e-6,'TolFun',1e-6); % Fminsearch choices
x2=fminsearchbnd(@(x) WILCOUmodelFull2(x,spec1,omegaN),[0.1 0.5 0.5 1],[0 0 0.01 0],[inf pi 1 inf],options); % Whittle likelihood eq (17) all frequencies
%% storing values from Whittle likelihood eq (17) all frequencies
alpha = x2(1); beta =  x2(2); rho = x2(3); A = x2(4);
[dum1,PF] = max(spec1); PF2 = NN+1-PF;
psi = 0.5*(angle(JZ(PF))+angle(JZ(PF2)));
psi = mod(psi+pi/2,pi)-pi/2;
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
PP(ii,:) = [alpha1e beta1e alpha2e beta2e A1e];
%% Bootstrap with periodogram
BB = 100; specS = spec1; BPS = zeros(BB,5);
for jj = 1:BB
specb = specS.*exprnd(1,1,NN);
x2B=fminsearchbnd(@(x) WILCOUmodelFull2(x,specb,omegaN),[0.1 0.5 0.5 1],[0 0 0.01 0],[inf pi 1 inf],options);
alpha = x2B(1); beta =  x2B(2); rho = x2B(3); A = x2B(4);
[dum1,PF] = max(spec1); PF2 = NN+1-PF;
psi = 0.5*(angle(JZ(PF))+angle(JZ(PF2)));
psi = mod(psi+pi/2,pi)-pi/2;
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
BPS(jj,:) = [alpha1e beta1e alpha2e beta2e A1e];
end
BTS1(ii,:) = std(BPS);
%% Bootstrap with Epanecnikov kernel
BB = 100; BPS = zeros(BB,5);
epan = 0.75*(1-(-1:0.1:1).^2); NE = length(epan);
specepan = zeros(1,NN);
for kk = 1:NN
    specepan(kk) = sum(spec1(mod([kk-(NE-1)/2:kk+(NE-1)/2]-1,NN)+1).*epan/sum(epan));
end
for jj = 1:BB
specb = specepan.*exprnd(1,1,NN);
x2B=fminsearchbnd(@(x) WILCOUmodelFull2(x,specb,omegaN),[0.1 0.5 0.5 1],[0 0 0.01 0],[inf pi 1 inf],options);
alpha = x2B(1); beta =  x2B(2); rho = x2B(3); A = x2B(4);
[dum1,PF] = max(spec1); PF2 = NN+1-PF;
psi = 0.5*(angle(JZ(PF))+angle(JZ(PF2)));
psi = mod(psi+pi/2,pi)-pi/2;
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
BPS(jj,:) = [alpha1e beta1e alpha2e beta2e A1e];
end
BTS2(ii,:) = std(BPS);
end
%% VALUS FOR TABLE 3
true = [alpha1 beta1 alpha2 beta2 A1];
100*std(PP)./abs(true(1:5)) % Monte Carlo
100*mean(BTS1)./abs(true(1:5)) % Bootstrap periodogram
100*mean(BTS2)./abs(true(1:5)) % Bootstrap Epanecnikov