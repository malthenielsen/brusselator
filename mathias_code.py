import numpy as np
from matplotlib import pyplot as plt
plt.style.use('science')
from scipy.integrate import RK45,solve_ivp
import numpy as np
#  from numba import njit, prange
import time
import matplotlib.pyplot as plt
#  from Functions_Simulations import *
from scipy.signal import find_peaks

#  def func():
#      tau0 = .1
#      b1 = 3
#
#      b3 = np.random.random()*5
#
#      b4 = np.random.random()*10
#
#      b2 = b3*b1**2/b4**2 - tau0
#      return b2
#
#  bin = []
#  for i in range(10000):
#      bin.append(func())
#
#  fig, ax = plt.subplots(1,1, figsize = (8,6))
#  ax.hist(bin, bins = np.linspace(-.5,3,500))
#  plt.show()
#  exit()




np.random.seed(42)
Amps = []
OM = []


for test2 in range(3,33):


    tau0 = .1

    b2 = 1000

    while b2>100:

        b1 = 3

        b3 = np.random.random()*5

        b4 = np.random.random()*10

        b2 = (b4 + b3*b1**2/b4**2 - tau0)


    om0 = 0.5*np.sqrt(4*b3*b1**2-tau0**2)


    omega = 0.86

    Tmax = 500

    acc = 0

    for test in range(1):

        A1 = 0.01

        b2 = (b4 + b3*b1**2/b4**2 - tau0)# + #0.01*test2

        mat = np.array([[b2 - b4, b3*(b1**2/b4**2)],[-b2 , -b3*(b1**2/b4**2)]])
        print(mat)
        print(np.linalg.eigvals(mat))



        #  b2 = b4

        #  print(test2,om0,b1,b2,b3,b4)
        print(b2,b4,b3*b1**2/b4**2 - tau0)
        print(b1, b3)

        ts = np.linspace(0, Tmax, Tmax*1000)


        bru1 =  lambda T,Y: [b1*(1+A1*np.sin(omega*T)) - (b2+b4)*Y[0] + b3*Y[0]**2*Y[1],

                             b2*Y[0] - b3*Y[0]**2*Y[1]]


        x0 = 0*b1/b4

        y0 = 0*b2*b4/(b3*b1)

        print(x0, y0)

        sol = solve_ivp (bru1, [0, Tmax], [x0, y0],t_eval=ts,rtol=1e-10)


        T = sol.t

        Y = sol.y

        t = T; y0 = Y[:][0]; y1 = Y[:][1];

        n_sig= 2

        x0 = b1/b4

        y00 = b2*b4/(b3*b1)


        xx = np.linspace(x0 - x0/20, x0 + x0/20, 20)
        yy = np.linspace(y00 - y00/20, y00 + y00/20, 20)
        #  xx = np.linspace(0,15,20)
        #  yy = np.linspace(0,20,20)

        #  for h in range(30):
        #
        #
        #      fig, ax = plt.subplots(1,1, figsize = (10,8))
        #      for k in range(20):
        #          for l in range(20):
        #              x0 = xx[k]
        #              y0 = yy[l]
        #              ts = np.linspace(0, Tmax, Tmax*3)
        #              #  bru1 =  lambda T,Y: [b1*(1+A1*np.sin(omega*T)) - (b2+b4)*Y[0] + b3*Y[0]**2*Y[1],
        #                                   #  b2*Y[0] - b3*Y[0]**2*Y[1]]
        #              #  print(x0, y0)
        #              #  sol = solve_ivp (bru1, [0, Tmax], [x0, y0],t_eval=ts,rtol=1e-3)
        #              #  Y = sol.y
        #              dx = b1*(1+A1*np.sin(omega*h*2)) - (b2+b4)*x0 + b3*x0**2*y0
        #              dy = b2*x0 - b3*x0**2*y0
        #              #  dx = Y[:][0][1000] - Y[:][0][999]
        #              #  dy = Y[:][1][1000] - Y[:][1][999]
        #
        #              ax.quiver(x0, y0,dx/10, dy/10)
        #
        #      ax.plot(Y[:][0], Y[:][1])
        #      x0 -= .5
        #      ax.set_xlim(min(xx), max(xx))
        #      ax.set_ylim(min(yy), max(yy))
        #      #  ax.set_xlim(-2,15)
        #      #  ax.set_ylim(-3,25)
        #      plt.savefig(f'fig_{h}')
        #  print('done')
        #      #  plt.show()
        #
        #

        
        #  peaks,valleys = Find_Peaks(y0,n_sig)
        #  peaks, _ = find_peaks(y0)
        #  valleys, _ = find_peaks(-y0)
        #
        #  Amp = 0
        #
        #  if (len(peaks)>3 and len(valleys)>3):
        #
        #      Top = y0[peaks[-3:]]
        #
        #      Bot = y0[valleys[-3:]]
        #
        #      Amp = np.mean(0.5*(Top-Bot))/np.mean(y0)
        #
        #  if (test == 0 and Amp<10000):
        #
        #      Amps.append(Amp)
        #
        #      OM.append(om0)
        #
        #  #  if (Amp > 1):
        #
        #      acc = 1;
        #
        #  #  fig, ax = plt.subplots(3,1, figsize = (10,6), sharex = True)
        #  fig, ax = plt.subplots(1,1, figsize = (10,6), sharex = True)
        #  #  ax[0].scatter(t,y0)
        #  ax.plot(y0,y1)
        #  #  ax[2].plot(t,A1*np.sin(omega*t))
        #  plt.show()
        #
        #
