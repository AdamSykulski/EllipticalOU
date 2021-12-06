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
PFULL = zeros(TT,5); PMARG = zeros(TT,5); PFULLSP = zeros(TT,5); PMARGSP = zeros(TT,5); % other initial entries
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
omegaN=0:2*pi/NN:2*pi*(1-1/NN); omegaN=fftshift(omegaN); omegaN(1:floor(NN/2))=omegaN(1:floor(NN/2))-2*pi; % Fourier frequencies
%%
JZ = (1/sqrt(NN))*fftshift(fft(ZNew4)); % Fourier transform
JZC = (1/sqrt(NN))*fftshift(fft(conj(ZNew4))); % Fourier transform of conjugate
spec1=(1/NN)*fftshift(abs(fft(ZNew4)).^2); % periodogram
spec2=(1/NN)*fftshift((fft(ZNew4)).*conj(fft(conj(ZNew4)))); % complementary spectrum
options=optimset('gradobj','on','MaxFunEval',10000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-10); % Fminsearch choices
x1=fminsearchbnd(@(x) WILCOUmodelFullF(x,JZ,JZC,omegaN,NN),[0.1 0.5 0.5 1 0],[0 0 0.01 0 -pi/2],[inf pi 1 inf pi/2],options); % Whittle likelihood eq (15) all frequencies
x2=fminsearchbnd(@(x) WILCOUmodelFull2(x,spec1,omegaN),[0.1 0.5 0.5 1],[0 0 0.01 0],[inf pi 1 inf],options); % Whittle likelihood eq (17) all frequencies
LF1 = sum(omegaN<betaP)-24; UF1 = sum(omegaN<betaP)+24; % 49 Fourier frequencies around the peak
LF2 = NN-UF1+1; % lower limit of frequencies used in semi-parametric fits
UF2 = NN-LF1+1; % upper limit of frequencies used in semi-parametric fits
x3=fminsearchbnd(@(x) WILCOUmodelRangeF(x,JZ,JZC,omegaN,LF1,UF1,NN),[0.1 0.5 0.5 1 0],[0 0 0.01 0 -pi/2],[inf pi 1 inf pi/2],options); % Whittle likelihood eq (15) reduced frequencies
x4=fminsearchbnd(@(x) WILCOUmodelRange2(x,spec1,omegaN,LF1,UF1,LF2,UF2),[0.1 0.5 0.5 1],[0 0 0.01 0],[inf pi 1 inf],options); % Whittle likelihood eq (17) reduced frequencies
%% storing values from Whittle likelihood eq (15) all frequencies
alpha = x1(1); beta =  x1(2); rho = x1(3); psi = x1(5); A = x1(4);
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
PFULL(ii,:) = [alpha1e beta1e alpha2e beta2e A1e];
%% storing values from Whittle likelihood eq (17) all frequencies
alpha = x2(1); beta =  x2(2); rho = x2(3); A = x2(4);
[dum1,PF] = max(spec1); PF2 = NN+1-PF;
psi = 0.5*(angle(JZ(PF))+angle(JZ(PF2))); % non-parametric estimate of psi
psi = mod(psi+pi/2,pi)-pi/2;
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
PMARG(ii,:) = [alpha1e beta1e alpha2e beta2e A1e];
%% storing values from Whittle likelihood eq (15) reduced frequencies
alpha = x3(1); beta =  x3(2); rho = x3(3); psi = x3(5); A = x3(4);
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
PFULLSP(ii,:) = [alpha1e beta1e alpha2e beta2e A1e];
%% storing values from Whittle likelihood eq (17) reduced frequencies
alpha = x4(1); beta =  x4(2); rho = x4(3); A = x4(4);
[dum1,PF] = max(spec1); PF2 = NN+1-PF;
psi = 0.5*(angle(JZ(PF))+angle(JZ(PF2))); % non-parametric estimate of psi
psi = mod(psi+pi/2,pi)-pi/2;
alpha1e = alpha;
beta1e = beta/2*(1/rho^2+rho^2);
alpha2e = beta/2*(rho^2-1/rho^2)*sin(2*psi);
beta2e = beta/2*(rho^2-1/rho^2)*cos(2*psi);
A1e = A*(1/rho^2+rho^2)/2;
PMARGSP(ii,:) = [alpha1e beta1e alpha2e beta2e A1e];
end
%% values for Table 2
true = [alpha1 beta1 alpha2 beta2 A1]
BIAS1 = 100*((mean(PFULL)-true)./abs(true))
BIAS2 = 100*((mean(PMARG)-true)./abs(true))
BIAS3 = 100*((mean(PFULLSP)-true)./abs(true))
BIAS4 = 100*((mean(PMARGSP)-true)./abs(true))
RMSE1 = 100*(sqrt(mean((PFULL-true).^2))./abs(true))
RMSE2 = 100*(sqrt(mean((PMARG-true).^2))./abs(true))
RMSE3 = 100*(sqrt(mean((PFULLSP-true).^2))./abs(true))
RMSE4 = 100*(sqrt(mean((PMARGSP-true).^2))./abs(true))
%% Generation of Figure 2
MC1 = figure; kde(100*((PMARG(:,1)-true(1))./true(1)))
hold on;  kde(100*((PMARG(:,2)-true(2))./true(2)))
hold on;  kde(100*((PMARG(:,3)-true(3))./abs(true(3))))
hold on;  kde(100*((PMARG(:,4)-true(4))./abs(true(4))))
hold on;  kde(100*((PMARG(:,5)-true(5))./true(5)))
legend('\alpha_1','\beta_1','\alpha_2','\beta_2','\sigma^2'); xlim([-50 50])
set(gca,'fontsize', 13); xlabel('percentage of true value')
saveas(MC1,'WILCOU_MC1.eps','epsc')