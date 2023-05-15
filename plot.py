import numpy as np
from matplotlib import pyplot as plt
plt.style.use('K_PAPER')

data = np.genfromtxt('scatter.txt' , delimiter = ',')
print(data)

fig, ax = plt.subplots(1,1, figsize = (6,6))
ax.scatter(data[:,0] + .1, data[:,1], alpha = .1)
ax.set_xlim(-2, 20)
ax.set_ylim(-10, 100)
ax.set_xlabel(r'$\frac{c a^2}{d^2}$')
ax.set_ylabel('Amp')
plt.show()


#  data = np.genfromtxt('results.txt' , delimiter = ',')
#
#  fig, ax = plt.subplots(2,1, figsize = (8,6))
#  ax[1].plot(data[:,1])
#  ax[0].plot(data[:,0])
#
#  fig, ax = plt.subplots(1,1, figsize = (6,6))
#  ax.plot(data[:,0], data[:,1])
#  #  plt.plot(data[:,1])
#  plt.show()

