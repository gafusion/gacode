from pygacode import expro
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,5)) ; ax = fig.add_subplot(111)
ax.set_xlabel('rho')
ax.set_ylabel('Bunit')

expro.expro_read('input.gacode')

x = expro.expro_rho
y = expro.expro_bunit

ax.plot(x,y)

expro.expro_read('input.gacode.g')

x = expro.expro_rho
y = expro.expro_bunit

ax.plot(x,y)

plt.show()

