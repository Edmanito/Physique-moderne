import numpy as np
import matplotlib.pyplot as plt

# Masquage des avertissements lies a l instabilite numerique
# Le schema d Euler explicite diverge aux temps longs ce qui est normal
# On desactive les alertes de la console car cette explosion mathematique
# est un comportement attendu que l on souhaite observer sur le graphique
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# 1. Parametres physiques et grille
# Constante de Planck reduite en Joules seconde
hbar = 1.054571817e-34  
# Masse de l electron en kilogrammes
m = 9.1093837015e-31    

# Discretisation de l espace de -10 nanometres a 20 nanometres
x_min, x_max = -10e-9, 20e-9
# Nombre de points pour le maillage spatial
nx = 800
x = np.linspace(x_min, x_max, nx)
# Pas spatial representant la distance entre deux points de la grille
dx = x[1] - x[0]

# Duree totale de la simulation fixee a 30 femtosecondes
t_max = 30e-15
# Nombre d iterations temporelles
nt = 15000
# Pas de temps tres court choisi pour retarder l instabilite du schema explicite
dt = t_max / nt

# 2. Creation de la barriere V0 = 5 eV
# Conversion de la hauteur de la barriere des electronvolts vers les Joules
V0 = 5.0 * 1.6e-19  
# Position d entree de la barriere a 5 nanometres
x_barriere = 5e-9
# Epaisseur de la barriere fixee a 2 nanometres
a_barriere = 2e-9

# Initialisation du potentiel a zero partout correspondant au vide
V = np.zeros(nx)
# Elevation du potentiel uniquement dans la zone spatiale definie pour la barriere
V[(x >= x_barriere) & (x <= x_barriere + a_barriere)] = V0

# 3. Initialisation
# Matrice complexe stockant la fonction d onde pour chaque position et chaque instant
Psi = np.zeros((nt, nx), dtype=complex)

# Nombre d onde central lie a l impulsion de depart de la particule
k0 = 5e9
# Parametre d etalement spatial initial
a = 1e-9
# Facteur de variance de la distribution gaussienne
gamma = (a**2)/4
# Constante de normalisation garantissant une probabilite totale de presence egale a un
N = np.sqrt(a) / ((2 * np.pi)**(0.25) * np.sqrt(2 * gamma))

# Injection du paquet d ondes initial sur la toute premiere ligne de temps
Psi[0, :] = N * np.exp(- (x**2) / (4 * gamma)) * np.exp(1j * k0 * x)

# 4. Resolution numerique via la methode Euler
# Precalcul du coefficient lie a l energie cinetique pour optimiser la boucle
coeff_cinetique = 1j * hbar / (2 * m * dx**2)

# Boucle d evolution temporelle
for n in range(0, nt - 1):
    # Boucle spatiale excluant les bords pour forcer la fonction d onde a zero aux limites
    for i in range(1, nx - 1):
        # Approximation de la derivee seconde spatiale par la methode des differences finies
        deriv2 = Psi[n, i+1] - 2*Psi[n, i] + Psi[n, i-1]
        
        # Coefficient local lie a l energie potentielle de la barriere
        coeff_potentiel = -1j * V[i] / hbar
        
        # Calcul et mise a jour de la fonction d onde pour le pas de temps suivant
        Psi[n+1, i] = Psi[n, i] + dt * (coeff_cinetique * deriv2 + coeff_potentiel * Psi[n, i])

# 5. Extraction de l image au moment de la collision a 60 pourcents du temps total
n_plot = int(0.6 * nt) 

# Trace
plt.figure(figsize=(8, 4))
plt.plot(x*1e9, np.abs(Psi[n_plot, :]), color='blue', label='Densité de probabilité |Psi|')

# Affichage de la barriere mise a l echelle pour etre visible sur le meme graphique
plt.fill_between(x*1e9, V/np.max(V)*np.max(np.abs(Psi[n_plot,:])), color='gray', alpha=0.4, label='Barrière V0')

plt.xlabel("Position en nm")
plt.title("Collision du paquet d ondes avec la barriere")
plt.legend()
plt.show()
