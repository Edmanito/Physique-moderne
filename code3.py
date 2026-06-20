import numpy as np
import matplotlib.pyplot as plt

# 1. Parametres de la barriere
# Tableau contenant les differentes epaisseurs de barriere a etudier allant de 0.1 a 5 nanometres
largeurs = np.linspace(0.1e-9, 5e-9, 100)

# 2. Temps de traversee classique dit temps libre
# Vitesse de groupe du paquet d ondes estimee en metres par seconde
vitesse_groupe = 5.8e5 
# En physique classique le temps de vol est directement proportionnel a la distance a parcourir
temps_libre = largeurs / vitesse_groupe

# 3. Temps de traversee quantique et phenomene de saturation
# Valeur asymptotique maximale du temps de traversee en secondes
tau_max = 2.5e-15 
# Constante d attenuation de l onde evanescente a l interieur de la barriere
kappa = 2e9 
# Modele analytique montrant que le temps de traversee stagne pour les barrieres larges
temps_tunnel = tau_max * (1 - np.exp(-2 * kappa * largeurs))

# 4. Visualisation graphique comparant les deux comportements
plt.figure(figsize=(8, 4))

# Conversion des largeurs en nanometres et des temps en femtosecondes pour la lisibilite du graphique
plt.plot(largeurs*1e9, temps_libre*1e15, 'r--', label='Temps libre classique')
plt.plot(largeurs*1e9, temps_tunnel*1e15, 'b-', linewidth=2, label='Temps tunnel quantique')

plt.xlabel("Largeur de la barriere a (nm)")
plt.ylabel("Temps de traversee en fs")
plt.title("Phenomene de saturation du temps")
plt.grid(True, alpha=0.3)
plt.legend()

plt.show()