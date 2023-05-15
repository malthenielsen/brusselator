import numpy as np
from matplotlib import pyplot as plt
plt.style.use('science')
np.random.seed(42)
from scipy.integrate import RK45,solve_ivp

Amps = []

OM = []

Pars = []

Tmax = 100

ts = 0.01

for test2 in range(28,1000):

    print(test2)

    A1 = 0.1

    omega = 1.0

    lamR = -.1

    acc = 0

    while acc == 0:

        b2 = A1 + np.random.random()

        b3 = 4*np.random.random()

        b4 = 10*np.random.random()


        beta = 2*lamR + b3 + 2*b2**2*b4/b3**2

        gamma = 2*lamR*b2**2*b4/b3**2 - b2**2*b4/b3 + b4**2*b2**4/b3**4

        Det = beta**2-4*gamma

        if (Det>0):

            b1 = 0.5*(-beta + np.sqrt(Det))

            if (b1>0):

                acc = 1


    ts = np.linspace(0, Tmax, Tmax*300)

    bru1 =  lambda T,Y: [-b3*Y[0] + b1*Y[1] + b4*Y[0]**2*Y[1],

                          A1*np.sin(omega*T) + b2 - b1*Y[1] - b4*Y[0]**2*Y[1]]


    x0 = b2/b3

    y0 = b2*b3**2/(b1*b3**2 + b2**2*b4)

    sol = solve_ivp (bru1, [0, Tmax], [x0+1, y0+1],t_eval=ts,rtol=1e-5)
    T = sol.t

    Y = sol.y

    t = T; y0 = Y[:][0]; y1 = Y[:][1];
    fig, ax = plt.subplots(3,1, figsize = (10,6), sharex = True)
    ax[0].plot(t,y0)
    ax[1].plot(t,y1)
    ax[2].plot(t,A1*np.sin(omega*t))
    plt.show()



