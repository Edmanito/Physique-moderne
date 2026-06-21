import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques
hbar = 1.054571817e-34
m = 9.1093837015e-31
eV = 1.602176634e-19

# Parametres de la barriere
V0_eV = 5.0
V0 = V0_eV * eV
a = 2e-9

# Creation d un tableau d energies allant de tres faible a 10 electronvolts
energies_eV = np.linspace(0.1, 10, 500)
energies = energies_eV * eV

# Initialisation du tableau pour stocker les probabilites de transmission
T = np.zeros(len(energies))

# Calcul analytique issu de l etude des etats stationnaires
for i, E in enumerate(energies):
    if E < V0:
        # Regime tunnel avec attenuation exponentielle sous la barriere
        kappa = np.sqrt(2 * m * (V0 - E)) / hbar
        # Formule exacte du coefficient de transmission pour ce regime
        denominateur = 1 + (V0**2 * np.sinh(kappa * a)**2) / (4 * E * (V0 - E))
        T[i] = 1 / denominateur
    else:
        # Regime oscillatoire au dessus de la barriere
        k_prime = np.sqrt(2 * m * (E - V0)) / hbar
        # Formule exacte utilisant les ondes propagatives
        denominateur = 1 + (V0**2 * np.sin(k_prime * a)**2) / (4 * E * (E - V0))
        T[i] = 1 / denominateur

# Trace du graphique
plt.figure(figsize=(8, 4))

# Trace principal de la transmission
plt.plot(energies_eV, T, color='purple', linewidth=2, label='Transmission quantique')

# Ajout d une ligne verticale pour marquer la limite classique de la barriere
plt.axvline(x=V0_eV, color='red', linestyle='--', label='Hauteur de la barriere V0')

plt.xlabel("Energie de la particule E en eV")
plt.ylabel("Coefficient de transmission T")
plt.title("Probabilite de transmission en fonction de l energie")
plt.grid(True, alpha=0.3)
plt.legend()

plt.show()
