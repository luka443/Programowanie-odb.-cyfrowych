import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Initialize parameters
from main import Nx

inputType = 2
Kiter = 10
N = 16
P = 4
h = [ 1, 0.5, 0, -0.25 ]
Nh = len(h)
step = (Nx - 1) / Nx  # Wrong sampling ratio

# Calculate channel frequency response
H = np.fft.fft( np.concatenate((h, np.zeros(N-Nh))) )/N

# Signal generation and transmission
x = []
S1ref = []

for k in range(Kiter):
    S0 = (2 * np.round(np.random.rand(int(N/2-1))) - 1) + 1j * (2 * np.round(np.random.rand(int(N/2-1))) - 1)
    S1 = np.concatenate(([0], S0, [0], np.conj(S0[::-1])))
    x1 = np.fft.ifft(S1)
    S1ref.extend(S1)
    x1p = np.concatenate((x1[N-P:], x1))
    x.extend(x1p.tolist())

y = np.convolve(x ,h)
y = y[0:Nx]

# Introduce wrong sampling ratio
y1 = interp1d(range(Nx), y, kind='cubic')
y = y1(np.arange(0, Nx - 1, step))

# Receiver
y1p = np.reshape(y, (N+P, Kiter))
y1 = y1p[P:N+P]
Y1e = np.fft.fft(y1) / N / np.reshape(H, (N, 1))

# Compute error
error2 = np.max(np.abs(np.array(S1ref) - Y1e.flatten()))
Y1e = Y1e[1:int(N/2)].tolist() + Y1e[int(N/2)+2:].tolist()
Y1e = np.array(Y1e).flatten()

# Plot the constellation diagram
plt.scatter(np.real(Y1e), np.imag(Y1e), color='b')
plt.grid()
plt.title('I=f(Q)')
plt.show()