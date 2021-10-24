import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import cmath
import math

omega = 1
delta = 2


def DifUr(variab_x_y, Axis_t):
    dxdt = variab_x_y[1]
    pow(omega, 2)
    dydt = -(2 * delta) * variab_x_y[1] - (math.pow(omega, 2)) * variab_x_y[0]
    return [dxdt, dydt]


Axis_t = np.linspace(0, 5)
variab_x_y_0 = [-1, 19]
variab_x_y = odeint(DifUr, variab_x_y_0, Axis_t)

dxdt = variab_x_y[:, 0]
dydt = variab_x_y[:, 1]
plt.plot(Axis_t, dxdt, '-o', Axis_t, dydt, '-o', linewidth=2.0)
plt.xlabel('t')
plt.ylabel('x(t)    y(t)')
plt.grid(True)
plt.show()

dydt = np.linspace(-2.0, 8.0, 20)
dxdt = np.linspace(-2.0, 2.0, 20)

Y1, Y2 = np.meshgrid(dydt, dxdt)

Axis_t = 0

u, v = np.zeros(Y1.shape), np.zeros(Y2.shape)

NI, NJ = Y1.shape

for i in range(NI):
    for j in range(NJ):
        x = Y1[i, j]
        y = Y2[i, j]
        yprime = DifUr([x, y], Axis_t)
        u[i, j] = yprime[0]
        v[i, j] = yprime[1]
Q = plt.quiver(Y1, Y2, u, v, color='r')
plt.xlabel('$y_1$')
plt.ylabel('$y_2$')
for y20 in [0, 0.5, 1, 1.5, 2, 2.5]:
    tspan = np.linspace(0, 50, 200)
    y0 = [0.0, y20]
    ys = odeint(DifUr, y0, tspan)
    plt.plot(ys[:, 0], ys[:, 1], 'b-')
    plt.plot([ys[0, 0]], [ys[0, 1]], 'o')
    plt.plot([ys[-1, 0]], [ys[-1, 1]], 's')
plt.xlim([-2, 8])
plt.ylim([-4, 4])

plt.show()

# ------------------------------------------------------------------------------------------#

# ! Вычисление АФЧХ, АЧХ, ФЧХ и их графики
# * Параметры системы
a = 2
i = 0
v = np.linspace(0, 5)

AFCH = (a) / ((omega ** 2) - (v ** 2) + (2 * delta * v * 1j))  # формула АФЧХ
ACH = np.absolute(AFCH)  #  По определию АЧХ это модуль от АФЧХ
FCH = []  # ФЧХ это аргумент от АФЧХ

# ! Заполнение массива значениями для ФЧХ
while i < len(v):
    fi = cmath.phase(AFCH[i])
    FCH.append(fi)
    i = i + 1

plt.plot(AFCH.real, AFCH.imag, 'o-', label='АФЧХ')
plt.xlabel('Re')
plt.ylabel('Im')
plt.legend()
plt.grid()
plt.show()

plt.plot(v, ACH, 'o-', label='АЧХ')
plt.xlabel('v')
plt.ylabel('A(v)')
plt.legend()
plt.grid()
plt.show()

plt.plot(v, FCH, 'o-', label='ФЧХ')
plt.xlabel('v')
plt.ylabel('fi(v)')
plt.legend()
plt.grid()
plt.show()
