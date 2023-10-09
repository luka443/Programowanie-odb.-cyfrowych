% lab18_ex_fm.m
% Example of FM modulation and demodulation
clear all; close all;

% Modulating signal m(t)
[x,fs] = audioread('GOODBYE.WAV'); % read speech from file
Nx = length(x);                    % number of signal samples
x = x(:,1)';                       % take only one channel, 1 or 2
dt=1/fs; t=dt*(0:Nx-1);            % time
df=1/(Nx*dt); f=df*(0:Nx-1);       % frequency
x = cos(2*pi*2*t);                % alternative modulating signal
%x = sin(2*pi*2*t)+ 0.5*sin(pi*t);
figure; plot(t,x); xlabel('t [s]'); grid; title('x(t)'); pause
figure; spectrogram(x,256,192,512,fs,'yaxis'); title('STFF of x(t)'); pause

% FM modulation
fc = 4000;                            % carrier frequency: 0 or 4000 Hz
BW = 0.9*fs;                       % available bandwidth
fmax = 3500;                       % maximum modulating frequency
df = (BW/(2*fmax)-1)*fmax,         % calculated frequency modulation depth
df = 1000;                         % arbitrary chosen freq modulation depth 
y = exp( j *2*pi*(fc*t + df*cumsum(x)*dt) );    % signal modulated in frequency
Y = abs(fft(y)/Nx);                             % its DFT spectrum
figure; plot(f,Y); grid; xlabel('f [Hz]'); title('|Y(f)|'); pause
figure; plot(f,20*log10(Y)); grid; xlabel('f [Hz]'); title('|Y(f)| [dB]'); pause
figure; spectrogram(y,256,192,512,fs,'yaxis); title('STFT of y(t)'); pause

% FM demodulation methods 
ang = unwrap(angle(y)); fi1 = 1/(2*pi)*(ang(2:end)-ang(1:end-1)) / dt;      % M1
fi2 = (1/(2*pi))*angle( y(2:Nx).*conj( y(1:Nx-1) ) ) / dt;                  % M2
fi3 = (1/(2*pi))*angle( y(3:Nx).*conj( y(1:Nx-2) ) ) / (2*dt); fi3=[fi3 0]; % M3
fi4 = (1/(2*pi))*...                                                        % M4
      (real(y(2:end-1)).*(imag(y(3:end))-imag(y(1:end-2)))-...              % M4
       imag(y(2:end-1)).*(real(y(3:end))-real(y(1:end-2))) )/(2*dt); fi4=[fi4 0];
fi5 = 1/(2*pi)*(real(y(1:end-1)).*imag(y(2:end))-imag(y(1:end-1)).*real(y(2:end)))/dt;
figure; nn = 1 : length(fi1);
plot(nn,fi1,'r',nn,fi2,'g',nn,fi3,'b',nn,fi4,'k',nn,fi5,'m'); title('Calculated angles'); 
xlabel('t [s]'); legend('\phi_1','\phi_2','\phi_3','\phi_4','\phi_5'); pause

xest = ( fi2 - fc ) / df ;        % recovered modulating signal
xest = xest(1:end-1);
x = x(2:end-1);
figure; plot(t(2:Nx-1),x,'r-',t(2:Nx-1),xest,'b-'); xlabel('t [s]');
title('Original and demodulated signal'); grid; legend('Origin','Demod'); pause
figure; spectrogram(xest,256,192,512,fs,'yaxis'); title('STFT of xest(t)'); pause

% ERROR after frequency MOD & DEMOD  
ERROR_SIGNAL = max( abs( x - xest ) ), pause  % FM demodulation error
soundsc(x,fs); pause                          % playing the original signal
soundsc(xest,fs); pause                       % playing the demodulated signal
