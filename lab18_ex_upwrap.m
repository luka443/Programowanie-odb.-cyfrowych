% lab18_ex_unwrap.m
% Example of signal angle wrapping and un-wrapping
clear all; close all;

n = 0 : 400;                                   % sample indexes
phi = 3*pi*sin(2*pi/200*n);                    % changing angle in radians
x = exp( j*phi );                              % generated complex-value signal:
                                               % cos(phi)+j*sin(phi)
phi1 = atan2( imag(x), real(x) );              % restoring the angle #1
phi2 = angle( x );                             % restoring the angle #2
error = max(abs( phi1 - phi2 )),               % error, should be equal 0
figure; plot(n,phi,'k--',n,phi1,'b'); grid;    % figure
xlabel('n'); ylabel('[rad]'); title('\phi (n): original and calculated'); pause

phi1 = unwrap( phi1 );                         % phase unwraping
figure; plot(n,phi,'k--',n,phi1,'b'); grid;    % figure
xlabel('n'); ylabel('[rad]'); title('\phi (n): original and calculated'); pause
