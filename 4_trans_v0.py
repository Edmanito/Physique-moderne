import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques fondamentales pour la mecanique quantique
hbar = 1.054571817e-34
m = 9.1093837015e-31
eV = 1.602176634e-19

# Parametres initiaux de la particule et epaisseur du mur
k0 = 5e9
a = 2e-9

# Calcul de l energie cinetique de la particule incidente
E = hbar**2 * k0**2 / (2*m)

# Creation d un tableau contenant des hauteurs de barriere de 1 a 6 electronvolts
V0_eV = np.linspace(1, 6, 100)

# Initialisation de la liste pour stocker les probabilites de transmission
T = []

# Boucle pour calculer la transmission pour chaque hauteur de barriere testee
for V in V0_eV:
    
    # Conversion de l energie de la barriere vers les Joules pour le calcul
    Vj = V * eV
    
    # Calcul du coefficient d attenuation de l onde evanescente a l interieur du mur
    kappa = np.sqrt(2*m*(Vj-E))/hbar
    
    # Application de l approximation classique de la barriere epaisse
    T.append(np.exp(-2*kappa*a))

# Configuration et affichage du graphique final
plt.figure(figsize=(8,4))

# Trace de la courbe decroissante representant l attenuation
plt.plot(V0_eV, T, color='green', linewidth=2)

# Etiquettes des axes clarifiees pour etre parfaitement lisibles
plt.xlabel("Hauteur de la barriere V0 en eV")
plt.ylabel("Coefficient de transmission T")
plt.title("Influence de la hauteur V0 sur la probabilite de transmission")
plt.grid(True, alpha=0.3)

plt.show()
