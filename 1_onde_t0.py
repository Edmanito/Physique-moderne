import numpy as np
import matplotlib.pyplot as plt
#Push de test
# 1. Paramètres physiques et grille spatiale

# Définition de l'axe spatial x (de -10 nm à 20 nm) avec 1000 points pour la résolution
x = np.linspace(-10e-9, 20e-9, 1000)
# Paramètres caractéristiques du paquet d'ondes initial
a = 1e-9         # Paramètre lié à l'étalement spatial (largeur) du paquet (en mètres)
k0 = 5e9         # Nombre d'onde central (lié à l'impulsion de la particule p = hbar * k0)
gamma = (a**2)/4 # Facteur de variance spatiale de la gaussienne

# 2. Construction mathématique de la fonction d'onde à t=0

# Calcul de la constante de normalisation N
# Elle garantit que la probabilité totale de présence vaut 1 (intégrale de |Psi|^2 dx = 1)
N = np.sqrt(a) / ((2 * np.pi)**(0.25) * np.sqrt(2 * gamma))

# Expression analytique du paquet d'ondes gaussien :
# Le paquet est le produit d'une enveloppe réelle  et d'une phase complexe
psi_0 = N * np.exp(-(x**2)/(4*gamma)) * np.exp(1j * k0 * x)

# 3. Tracé du graphique

plt.figure(figsize=(8, 4))

# Affichage de la partie réelle 
# On multiplie x par 1e9 pour convertir les mètres en nanomètres sur l'axe des abscisses
plt.plot(x*1e9, np.real(psi_0), label='Partie réelle', color='blue', alpha=0.6)

# Affichage de l'enveloppe globale (module de la fonction d'onde quantique)
plt.plot(x*1e9, np.abs(psi_0), 'r-', linewidth=2, label='Enveloppe |Psi|')

plt.xlabel("Position (nm)")
plt.ylabel("Amplitude")
plt.title("Paquet d'ondes gaussien initial (t=0)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()
