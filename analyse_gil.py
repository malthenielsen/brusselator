import numpy as np
from matplotlib import pyplot as plt
plt.style.use('science')


data = np.fromfile('results.bin', dtype = np.float64).reshape(3,-1)
print(data.shape)
peaks = np.fromfile('peaks.bin', dtype = np.int32)

fig,ax = plt.subplots(1,1, figsize = (8,6))
ax.plot(data[-1,:], data[0,:], color = 'black')
ax.scatter(data[-1,peaks], data[0,peaks], color = 'red')
#  ax.plot(data[-1,:], data[1,:])
#  ax.plot(data[0,:], data[1,:])
plt.show()

