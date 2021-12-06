%% INPUTS
Delta = 1/1000; % Euler Maruyama increment
N = 1000; % Time Series Length
alpha1 = 0.02; beta1 = 1; alpha2 = -0.5; beta2 = -0.3; A1 = 2; % elliptical OU parameters (FEEL FREE TO TRY DIFFERENT VALUES)
%% GENERATION CODE
rng(1); % seed
R1 = A1*(-beta2-1i*alpha2)/beta1; % Noise relation (equation at bottom of page 4 in paper)
T = N/Delta; % total number of increments in Euler Maruyama generation
ZZ = zeros(1,T);
for t = 1:T-1 % Euler Maruyama generation
    s1 = (A1+real(R1))/2; s2 = (A1-real(R1))/2; % magnitude of noisy in x/y components
    corr=imag(R1)/(2*sqrt(s1)*sqrt(s2)); % correlation between x/y noise components (please refer to Sykulski and Percival (2016) IEEE paper for details)
    B1 = randn(1,1); C1 = randn(1,1); D1 = corr*B1 + sqrt(1-corr^2)*C1; % correlated noise in x/y
    N1 = sqrt(s1)*B1; N2 = sqrt(s2)*D1; % scaled to correct amplitudes
    ZZ(t+1) = ZZ(t) - (alpha1-1i*beta1)*ZZ(t)*Delta - (alpha2-1i*beta2)*conj(ZZ(t))*Delta + sqrt(Delta)*(N1+1i.*N2); % Euler Maruyama forwards equation
end
ZNew3=ZZ(1/Delta:1/Delta:end); % subsample to correct sampling frequency
SZC = (1/N)*fftshift(abs(fft(ZNew3)).^2); % peridogoram
%% Trajectory plot
F1a=figure; plot3(Delta*(900001:1000000),real(ZZ(end-99999:end)),imag(ZZ(end-99999:end)),'k'); grid('on');
ylim([-10 10]); xlim([900 1000]); zlim([-10 10]);
xlabel('t'); ylabel('x(t)'); zlabel('y(t)');
hold on; plot3(Delta*(900001:1000000),10*ones(1,100000),imag(ZZ(end-99999:end)),'color',[.7 .7 .7]);
hold on; plot3(Delta*(900001:1000000),real(ZZ(end-99999:end)),-10*ones(1,100000),'color',[.7 .7 .7]);
set(gca,'fontsize', 16)
saveas(F1a,'WILCOU_EM1a.eps','epsc')
%% Spectrum plot
alpha = alpha1; % parameters in bivariate setting using Table 1 of paper
rho = ((abs(beta1)-sqrt(alpha2^2+beta2^2))/(abs(beta1)+sqrt(alpha2^2+beta2^2)))^0.25;
beta = sign(beta1)*sqrt(beta1^2-beta2^2-alpha2^2);
psi = -sign(beta1)*0.5*atan2(alpha2,-sign(beta1)*beta2);
A = A1*sqrt(beta1^2-beta2^2-alpha2^2)/abs(beta1);
omegaN=0:2*pi/N:2*pi*(1-1/N); omegaN=fftshift(omegaN); omegaN(1:floor(N/2))=omegaN(1:floor(N/2))-2*pi; % Fourier frequencies
SZZE = ((1/rho+rho)^2/4)*A./(alpha^2+(omegaN-beta).^2) + ((1/rho-rho)^2/4)*A./(alpha^2+(omegaN+beta).^2); % Spectrum (unaliased)
F2a=figure; plot(omegaN,10*log10(SZC),'k'); % plot of periodogram
hold on; plot(omegaN,10*log10(SZZE),'r','linewidth',2); xlim([-pi pi]); xlabel('\omega'); grid on; % overlaid plot of spectrum (unaliased)
K=10; omegaM = zeros(2*K+1,length(omegaN)); % next few lines calculate the aliased spectral density
for kk = 1:K
omegaM(K+1-kk,:) = omegaN-2*pi*kk;
omegaM(K+1+kk,:) = omegaN+2*pi*kk;
end
omegaM(K+1,:) = omegaN;
SZZ = sum(((1/rho+rho)^2/4)*A./(alpha^2+(omegaM-beta).^2) + ((1/rho-rho)^2/4)*A./(alpha^2+(omegaM+beta).^2));
hold on; plot(omegaN,10*log10(SZZ),'g','linewidth',2); xlim([-pi pi]); xlabel('\omega'); grid on; ylim([-30 50]); % overlaid plot of aliased spectrum
set(gca,'fontsize', 16)
saveas(F2a,'WILCOU_EM2a.eps','epsc')