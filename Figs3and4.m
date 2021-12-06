%% Preliminaries
load EarthWobble.mat % load data
TT = (1845+(363/365)):0.1:(2021+10/12); % Sampled time points
NN = length(Z); % length
omega=0:2*pi/NN:2*pi*(1-1/NN); omega=fftshift(omega); omega(1:floor(NN/2))=omega(1:floor(NN/2))-2*pi; % Fourier freuquencies
spec1=(1/NN)*fftshift(abs(fft(Z)).^2); % periodogram
%% DATA (Figure 3 Left)
F1=figure; plot3(TT,real(Z),imag(Z),'k'); grid('on');
ylim([-600 600]); zlim([-600 600]); xlim([1845 2022]);
xlabel('year'); ylabel('x (mas)'); zlabel('y (mas)');
hold on; plot3(TT,600*ones(1,NN),imag(Z),'color',[.7 .7 .7]);
hold on; plot3(TT,real(Z),-600*ones(1,NN),'color',[.7 .7 .7]);
set(gca,'fontsize', 13)
saveas(F1,'WILCOU_Fig1.eps','epsc')
%% Periodogram - with regions highlighted (Figure 3 Right)
F2=figure; plot(5*omega/pi,10*log10(spec1),'k'); xlim([-5 5]); grid on
LF1 = 709;
UF1 = 757;
LF2 = NN-UF1+2;
UF2 = NN-LF1+2;
spec2 = spec1(LF1:UF1);
omega2 = omega(LF1:UF1);
spec3 = spec1(LF2:UF2);
omega3 = omega(LF2:UF2);
hold on; plot(5*omega2/pi,10*log10(spec2),'b')
hold on; plot(5*omega3/pi,10*log10(spec3),'b')
LF1 = 699;
UF1 = 709;
LF2 = NN-UF1+2;
UF2 = NN-LF1+2;
spec2 = spec1(LF1:UF1);
omega2 = omega(LF1:UF1);
spec3 = spec1(LF2:UF2);
omega3 = omega(LF2:UF2);
hold on; plot(5*omega2/pi,10*log10(spec2),'r')
hold on; plot(5*omega3/pi,10*log10(spec3),'r')
xlabel('cycles/year'); ylabel('dB'); ylim([0 80])
set(gca,'fontsize', 13)
saveas(F2,'WILCOU_Fig2.eps','epsc')
%% Filter Chandler Wobble (Figure 4 Left)
LF1 = 709;
UF1 = 757;
LF2 = NN-UF1+1;
UF2 = NN-LF1+1;
FFTC = fft(Z); FFTC = fftshift(FFTC);
FFTF = [zeros(LF1-1,1); FFTC(LF1:UF1); zeros(LF2-UF1-1,1); FFTC(LF2:UF2); zeros(NN-UF2,1)];
Z2 = ifft(fftshift(FFTF));
F3=figure; plot3(TT,real(Z2),imag(Z2),'b'); grid('on');
ylim([-300 300]); zlim([-300 300]); xlim([1844 2022]);
xlabel('year'); ylabel('x (mas)'); zlabel('y (mas)');
hold on; plot3(TT,300*ones(1,NN),imag(Z2),'color',[.7 .7 .7]);
hold on; plot3(TT,real(Z2),-300*ones(1,NN),'color',[.7 .7 .7]);
set(gca,'fontsize', 13)
saveas(F3,'WILCOU_Fig3.eps','epsc')
%% Filtered annual wobble (Figure 4 Right)
LF1 = 699;
UF1 = 709;
LF2 = NN-UF1+1;
UF2 = NN-LF1+1;
FFTC = fft(Z); FFTC = fftshift(FFTC);
FFTF = [zeros(LF1-1,1); FFTC(LF1:UF1); zeros(LF2-UF1-1,1); FFTC(LF2:UF2); zeros(NN-UF2,1)];
Z2 = ifft(fftshift(FFTF));
F4=figure; plot3(TT,real(Z2),imag(Z2),'r'); grid('on');
ylim([-150 150]); zlim([-150 150]); xlim([1844 2020]);
xlabel('year'); ylabel('x (mas)'); zlabel('y (mas)');
hold on; plot3(TT,150*ones(1,NN),imag(Z2),'color',[.7 .7 .7]);
hold on; plot3(TT,real(Z2),-150*ones(1,NN),'color',[.7 .7 .7]);
set(gca,'fontsize', 13)
saveas(F4,'WILCOU_Fig4.eps','epsc')