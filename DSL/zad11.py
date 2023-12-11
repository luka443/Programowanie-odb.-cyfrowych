import numpy as np
import matplotlib.pyplot as plt

# Setup parameters
from main import Nx

inputType = 2
Kiter = 10
N = 16
P = 4
h = [ 1, 0.5, 0, -0.25 ]
Nh = len(h)
narrow_band_interference_amplitude = 0.1  # New parameters for the narrow-band interference
narrow_band_interference_frequency = 11.5  # New parameters for the narrow-band interference

# Calculate channel frequency response
H = np.fft.fft( np.concatenate((h, np.zeros(N-Nh))) )/N

# Allocate arrays for the transmitted signals and carrier states
x = []
S1ref = []

for k in range(Kiter):
    S0 = (2 * np.round(np.random.rand(int(N/2-1))) - 1) + 1j * (2 * np.round(np.random.rand(int(N/2-1))) - 1)
    S1 = np.concatenate(([0], S0, [0], np.conj(S0[::-1])))
    x1 = np.fft.ifft(S1)
    S1ref.extend(S1)
    x1p = np.concatenate((x1[N-P:], x1))
    x.extend(x1p.tolist())

# Compute signal through the channel and add narrow-band interference
y = np.convolve(x ,h)
y = y[0:Nx] + narrow_band_interference_amplitude * np.sin(2 * np.pi / N * narrow_band_interference_frequency * np.arange(Nx))

# Receiver
y1p = np.reshape(y, (N+P, Kiter))
y1 = y1p[P:N+P]
Y1e = np.fft.fft(y1) / N / np.reshape(H, (N, 1))

# Adjust carrier states and compute error
S1ref = np.array(S1ref)
error2 = np.max(np.abs(S1ref - Y1e.flatten()))
Y1e = Y1e[1:int(N/2)].tolist() + Y1e[int(N/2)+2:].tolist()
Y1e = np.array(Y1e).flatten()

# Plot the constellation diagram
plt.scatter(np.real(Y1e), np.imag(Y1e), color='b')
plt.grid()
plt.title('I=f(Q)')
plt.show()
#The graph shows the received carrier constellations with a narrow-band interference added.
# The interference is introduced as a sinusoid which has its frequency exactly in the middle between the 10th and 11th carrier.
# As expected, the resulting constellation points are more dispersed due to the narrow-band interference.