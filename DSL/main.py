import numpy as np
import matplotlib.pyplot as plt

inputType = 2
Kiter = 10
N = 16; P = 4
nstd = 0.05
h = [ 1, 0.5, 0, -0.25 ]; Nh = len(h)

# channel frequency response
H = np.fft.fft( np.concatenate((h, np.zeros(N-Nh))) )/N

# Transmitter
Nx = Kiter*(N+P)
x = []; x1ref = []; S1ref = []

for k in range(Kiter):

    if(inputType==1): x1 = np.random.randn(N)
    else:
        S0=(2*np.round(np.random.rand(int(N/2-1)))-1) + 1j*(2*np.round(np.random.rand(int(N/2-1)))-1)
        S1  = np.concatenate(([0], S0, [0], np.conj( S0[::-1]) ))
        x1 = np.fft.ifft( S1 )
        S1ref.extend(S1)
    x1p = np.concatenate((x1[N-P:], x1))
    x.extend(x1p.tolist())
    x1ref.extend(x1.tolist())

# Signal through the channel
y = np.convolve(x ,h); y = y[0:Nx]
y = y + nstd*np.random.randn(Nx)

# Receiver
y1p = np.reshape(y, (N+P, Kiter))
y1 = y1p[ P:N+P ]
Y1e = np.fft.fft(y1)/N / np.reshape(H, (N, 1))

if(inputType==1):
   y1e = np.fft.ifft(Y1e)
   error1 = np.max( np.abs( x1ref - y1e.flatten() ))
else:
   error2 = np.max( np.abs( S1ref - Y1e.flatten() ))
   Y1e=Y1e[1:int(N/2)].tolist() + Y1e[int(N/2)+1:].tolist(); Y1e=np.array(Y1e).flatten()
   plt.scatter(np.real(Y1e[:]),np.imag(Y1e[:]), color='b')
   plt.grid()
   plt.title('I=f(Q)')
   plt.show()


   #The figure above presents a scatter plot that indicates the received carrier constellation when using the 4-QAM modulation scheme (Quadrature Amplitude Modulation). Each point on the graph corresponds to a carrier (a specific frequency) in the signal received at the end of the Digital Subscriber Line (DSL).

#Both the real and imaginary components of the carriers are shown on the x-axis and y-axis respectively giving us an idea of how the symbols (signal states) are being received. In an ideal scenario, these would be closely grouped around the four 4-QAM constellation points at (-1, -1), (-1, 1), (1, -1) and (1, 1). Any deviation from these constellation points indicates errors introduced due to noise, interference, or signal distortion.