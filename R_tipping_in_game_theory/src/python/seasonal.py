import numpy as np
from scipy.fft import fft, ifft, fftfreq, fftshift
import matplotlib.pyplot as plt

t = np.linspace(0.0, 10*np.pi, 1000)
f = 5.0
y = np.cos(2*np.pi*f*t)
#spectrum = fftshift(fft(np.cos(3*t)))
#w = fftshift(fftfreq(t.shape[-1]))
spectrum = fft(y)
w = fftfreq(t.shape[-1],d=(t[1]-t[0]))
y_hat = ifft(spectrum)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
ax1 = plt.subplot2grid((3,4), (0,0), colspan=4, rowspan=1, fig=fig)
ax1.plot(t, y)
ax2 = plt.subplot2grid((3,4), (1,0), colspan=4, rowspan=1, fig=fig)
ax2.plot(w, spectrum.real)
ax3 = plt.subplot2grid((3,4), (2,0), colspan=4, rowspan=1, fig=fig)
ax3.plot(t, y_hat)

plt.show()
