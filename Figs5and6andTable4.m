%% Preliminaries
load EarthWobble.mat % load data
TT = (1845+(363/365)):0.1:(2021+10/12); % Sampled time points
NN = length(Z); % length
omega=0:2*pi/NN:2*pi*(1-1/NN); omega=fftshift(omega); omega(1:floor(NN/2))=omega(1:floor(NN/2))-2*pi; % Fourier freuquencies
spec1=(1/NN)*fftshift(abs(fft(Z)).^2); % periodogram
spec2=(1/NN)*fftshift((fft(Z)).*conj(fft(conj(Z)))); % complementary spectrum
JZ = (1/sqrt(NN))*fftshift(fft(Z)); % Fourier transform
JZC = (1/sqrt(NN))*fftshift(fft(conj(Z))); % Fourier transform of conjugate
%% Figure 5 (left)
LF1 = 709; % semi-parametric window, lower Fourier frequency
UF1 = 757; % semi-parametric window, upper Fourier frequency
options=optimset('gradobj','on','MaxFunEval',10000,'TolX',1e-10,'TolFun',1e-10); % Fminsearch choices
QG = [0.08 -(2*pi)*365.25/(433*10) 100]; % starting values
x1a=fminsearchbnd(@(x) WILCOUmodelRangeP(x,spec1,omega',LF1,UF1),QG,[0 -pi 0],[inf 0 inf],options); %optimisation
Q = [x1a(1) x1a(2) x1a(3)]; % parameter values
SZZ = Q(3)./(Q(1)^2+(omega-Q(2)).^2); % model spectrum
Q(1)*10 % damping in units of per year
Q(2)*10/(2*pi) % peak frequency in units of cycles per year
Q(3) % amplitude parameter
F5a=figure; % figure with confidence bands
boundedline(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),[-10*log10(chi2inv(0.975,2)/2); 10*log10(chi2inv(0.025,2)/2)],'r')
plot(5*omega(LF1:UF1)/pi,10*log10(spec1(LF1:UF1)),'k','linewidth',2);
hold on; plot(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),'color','r','linewidth',2); xlim(5*[omega(LF1) omega(UF1)]/pi);
xlabel('cycles/year'); ylabel('dB'); grid on
set(gca,'fontsize', 13)
saveas(F5a,'WILCOU_Fig5a.eps','epsc')
%% Bootstrap confidence intervals for Figure 5 analysis
BB = 100; % set to 10000 in paper for one run
x1aB = zeros(BB,3);
for ii = 1:BB
specb = spec1.*exprnd(1,NN,1);
x1aB(ii,:)=fminsearchbnd(@(x) WILCOUmodelRangeP(x,specb,omega',LF1,UF1),QG,[0 -pi 0],[inf 0 inf],options);
end
[quantile(x1aB(:,1)*10,0.025) quantile(x1aB(:,1)*10,0.975)] % damping in units of per year
[quantile(x1aB(:,2)*10/(2*pi),0.025) quantile(x1aB(:,2)*10/(2*pi),0.975)] % peak frequency in units of cycles per year
[quantile(x1aB(:,3),0.025) quantile(x1aB(:,3),0.975)] % amplitude parameter
%% Figure 5 (right)
LF1 = 709; % semi-parametric window, lower Fourier frequency
UF1 = 757; % semi-parametric window, upper Fourier frequency
LF2 = NN-UF1+1; % semi-parametric window, lower Fourier frequency
UF2 = NN-LF1+1; % semi-parametric window, upper Fourier frequency
QG = [0.08 -(2*pi)*365.25/(433*10) 0.9 100]; % starting values
x1b=fminsearchbnd(@(x) WILCOUmodelRange2(x,spec1,omega',LF1,UF1,LF2,UF2),QG,[0 -pi 0 0],[inf 0 1 inf],options); % optimisation
Q = [x1b(1) x1b(2) x1b(3) x1b(4)]; % parameter values
SZZ = ((1/Q(3)+Q(3))^2/4)*Q(4)./(Q(1)^2+(omega-Q(2)).^2) + ((1/Q(3)-Q(3))^2/4)*Q(4)./(Q(1)^2+(omega+Q(2)).^2); % model spectrum
F5b=figure; % figure with confidence bands
subplot(1,2,1);
boundedline(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),[-10*log10(chi2inv(0.975,2)/2); 10*log10(chi2inv(0.025,2)/2)],'r')
plot(5*omega(LF1:UF1)/pi,10*log10(spec1(LF1:UF1)),'k','linewidth',2);
hold on; plot(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),'color','r','linewidth',2);
xlim(5*[omega(LF1) omega(UF1)]/pi);
xlabel('cycles/year'); ylabel('dB'); grid on; ylim([15 75])
set(gca,'fontsize', 13)
subplot(1,2,2);
boundedline(5*omega(LF2:UF2)/pi,10*log10(SZZ(LF2:UF2)),[-10*log10(chi2inv(0.975,2)/2); 10*log10(chi2inv(0.025,2)/2)],'r')
plot(5*omega(LF2:UF2)/pi,10*log10(spec1(LF2:UF2)),'k','linewidth',2);
hold on; plot(5*omega(LF2:UF2)/pi,10*log10(SZZ(LF2:UF2)),'color','r','linewidth',2);
xlim(5*[omega(LF2) omega(UF2)]/pi);
xlabel('cycles/year'); ylabel('dB'); grid on; ylim([15 75])
set(gca,'fontsize', 13)
saveas(F5b,'WILCOU_Fig5b.eps','epsc')
%% Figure 6
LF1 = 699; % semi-parametric window, lower Fourier frequency
UF1 = 709; % semi-parametric window, upper Fourier frequency
LF2 = NN-UF1+1; % semi-parametric window, lower Fourier frequency
UF2 = NN-LF1+1; % semi-parametric window, upper Fourier frequency
QG = [0.08 -(2*pi)/10 0.9 100]; % starting values
x1c=fminsearchbnd(@(x) WILCOUmodelRange2(x,spec1,omega',LF1,UF1,LF2,UF2),QG,[0 -(2*pi)/10 0 0],[inf -(2*pi)/10 1 inf],options); % optimisation
Q = [x1c(1) x1c(2) x1c(3) x1c(4)]; % parameter values
SZZ = ((1/Q(3)+Q(3))^2/4)*Q(4)./(Q(1)^2+(omega-Q(2)).^2) + ((1/Q(3)-Q(3))^2/4)*Q(4)./(Q(1)^2+(omega+Q(2)).^2); % model spectrum
F6=figure; % figure with confidence bands
subplot(1,2,1);
boundedline(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),[-10*log10(chi2inv(0.975,2)/2); 10*log10(chi2inv(0.025,2)/2)],'r')
plot(5*omega(LF1:UF1)/pi,10*log10(spec1(LF1:UF1)),'k','linewidth',2);
hold on; plot(5*omega(LF1:UF1)/pi,10*log10(SZZ(LF1:UF1)),'color','r','linewidth',2);
xlim(5*[omega(LF1) omega(UF1)]/pi);
xlabel('cycles/year'); ylabel('dB'); grid on; ylim([15 75])
set(gca,'fontsize', 13)
subplot(1,2,2);
boundedline(5*omega(LF2:UF2)/pi,10*log10(SZZ(LF2:UF2)),[-10*log10(chi2inv(0.975,2)/2); 10*log10(chi2inv(0.025,2)/2)],'r')
plot(5*omega(LF2:UF2)/pi,10*log10(spec1(LF2:UF2)),'k','linewidth',2);
hold on; plot(5*omega(LF2:UF2)/pi,10*log10(SZZ(LF2:UF2)),'color','r','linewidth',2);
xlim(5*[omega(LF2) omega(UF2)]/pi);
xlabel('cycles/year'); ylabel('dB'); grid on; ylim([15 75])
set(gca,'fontsize', 13)
saveas(F6,'WILCOU_Fig6.eps','epsc')
%% analysis of annual oscillation
pp1=JZ(704);
pp2=JZ(1056);
psi = 0.5*(angle(pp1)+angle(pp2)) % non-parametric orientation
2*sqrt(abs(pp1).*abs(pp2))./(abs(pp1)+abs(pp2)) % non-parametric eccentricity
alpha = Q(1); beta =  Q(2); rho = Q(3); A = Q(4); % parameter estimates
alpha1 = alpha*10; % converting to elliptical OU form of equation (1)
beta1 = beta/2*(1/rho^2+rho^2)*10/(2*pi);
alpha2 = beta/2*(rho^2-1/rho^2)*sin(2*psi)*10/(2*pi);
beta2 = beta/2*(rho^2-1/rho^2)*cos(2*psi)*10/(2*pi);
A1 = A*(1/rho^2+rho^2)/2;
R1 = A*(1/rho^2-rho^2)*exp(2*1i*psi)/2;
ecc = sqrt(1-rho^4);
%% Values for Table 4 with bootstrap confidence intervals
BB = 100; % number of bootstraps, set to 10000 in paper for one run
x1cB = zeros(BB,4);
for ii = 1:BB
specb = spec1.*exprnd(1,NN,1);
x1cB(ii,:)=fminsearchbnd(@(x) WILCOUmodelRange2(x,specb,omega',LF1,UF1,LF2,UF2),QG,[0 -(2*pi)/10 0 0],[inf -(2*pi)/10 1 inf],options);
end
alphae = x1cB(:,1); betae =  x1cB(:,2); rhoe = x1cB(:,3); Ae = x1cB(:,4); ecce = sqrt(1-rhoe.^4);
alpha1e = alphae*10;
beta1e = betae/2.*(1./rhoe.^2+rhoe.^2)*10/(2*pi);
alpha2e = betae/2.*(rhoe.^2-1./rhoe.^2)*sin(2*psi)*10/(2*pi);
beta2e = betae/2.*(rhoe.^2-1./rhoe.^2)*cos(2*psi)*10/(2*pi);
A1e = Ae.*(1./rhoe.^2+rhoe.^2)/2;
% OUTPUT FOR TABLE 4
[alpha1 quantile(alpha1e,0.025) quantile(alpha1e,0.975)]
[beta1 quantile(beta1e,0.025) quantile(beta1e,0.975)]
[alpha2 quantile(alpha2e,0.025) quantile(alpha2e,0.975)]
[beta2 quantile(beta2e,0.025) quantile(beta2e,0.975)]
[A1 quantile(A1e,0.025) quantile(A1e,0.975)]
[ecc quantile(ecce,0.025) quantile(ecce,0.975)]