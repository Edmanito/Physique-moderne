import numpy as np
import matplotlib.pyplot as plt

hbar = 1.054571817e-34
m = 9.1093837015e-31
eV = 1.602176634e-19

k0 = 5e9
a = 2e-9

E = hbar**2 * k0**2 / (2*m)

V0_eV = np.linspace(1, 6, 100)

T = []

for V in V0_eV:

    Vj = V * eV

    kappa = np.sqrt(2*m*(Vj-E))/hbar

    T.append(np.exp(-2*kappa*a))

plt.figure(figsize=(8,4))
plt.plot(V0_eV, T)

plt.xlabel("Hauteur de barrière V0 (eV)")
plt.ylabel("Transmission T")
plt.title("Influence de V0 sur la transmission")
plt.grid(True)

plt.show()