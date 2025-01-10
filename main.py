import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import product
from fractions import Fraction



"""
################
### PARTIE 1 ###
################
"""

# ETAPE 1 : discretisation 
# Discrétisation pour la population de proies (x)
#1) : formulation discrète de la dérivée
# x'[i] = dx/dt = (x[i+1] - x[i]) / step

# 2) : injection de la formulation discrète dans l'équation
# x'[i] = dx/dt = alpha * x[i] - beta * x[i] * y[i]      ( apres avoir developper)
# (x[i+1] - x[i]) / step = alpha * x[i] - beta * x[i] * y[i]

# 3) : isolation de x[i+1]
# x[i+1] = x[i] + step * (alpha * x[i] - beta * x[i] * y[i])  



# Formes discretes : 
# x[i+1] = x[i] + step * (alpha * x[i] - beta * x[i] * y[i])  Lapins
# y[i+1] = y[i] + step * (delta * x[i] * y[i] - gamma * y[i]) Renards

# Initialisation des variables
STEP = 0.001
DAYS = 1000

# Paramètres du modèle
ALPHA, BETA, DELTA, GAMMA = 1/3, 1/3, 1/3, 1/3

# Conditions initiales
lapin = np.zeros(DAYS)
renard = np.zeros(DAYS)
time = np.zeros(DAYS)
lapin[0] = 1000
renard[0] = 2000

# Algorithme d'Euler
for i in range(1, DAYS):
    lapin[i] = lapin[i - 1] + STEP * (ALPHA * lapin[i - 1] - BETA * lapin[i - 1] * renard[i - 1])
    renard[i] = renard[i - 1] + STEP * (DELTA * lapin[i - 1] * renard[i - 1] - GAMMA * renard[i - 1])
    time[i] = time[i - 1] + STEP

# Tracer les résultats
plt.figure(figsize=(15, 6))
plt.plot(time, lapin, label="Lapins", linewidth=2)
plt.plot(time, renard, label="Renards", linewidth=2)
plt.xlabel("Temps")
plt.ylabel("Population")
plt.title("Modèle proie-prédateur de Lotka-Volterra (Lapins et Renards)")
plt.legend()
plt.grid()
plt.show()

"""
################
### PARTIE 2 ###
################
"""
# 1)
# Charger les données
DATA_FILE = 'populations_lapins_renards.csv'
df = pd.read_csv(DATA_FILE)

# Initialisation des données
days_data = df.shape[0]
lapin_reel = df['lapin'].values.astype(float)
renard_reel = df['renard'].values.astype(float)

# 2)
# Fonction pour simuler le modèle avec la methode d'euler avec des paramètres donnés
def simulate_model_euler(alpha, beta, delta, gamma, lapin0, renard0, step, days):
    lapin_theory = np.zeros(days)
    renard_theory = np.zeros(days)
    time = np.zeros(days)

    lapin_theory[0] = lapin0 / 1000
    renard_theory[0] = renard0 / 1000

    for i in range(1, days):
        lapin_theory[i] = lapin_theory[i - 1] + step * (
            alpha * lapin_theory[i - 1] - beta * lapin_theory[i - 1] * renard_theory[i - 1]
        )
        renard_theory[i] = renard_theory[i - 1] + step * (
            delta * lapin_theory[i - 1] * renard_theory[i - 1] - gamma * renard_theory[i - 1]
        )
        time[i] = time[i - 1] + step
    print((lapin_theory * 1000).dtype)
    return lapin_theory * 1000, renard_theory * 1000, time

# Fonction pour calculer l'erreur quadratique moyenne (MSE)
def calculate_mse(alpha, beta, delta, gamma, lapin_reel0, renard_reel0):
    lapin_pred, renard_pred, _ = simulate_model(
        alpha, beta, delta, gamma, lapin_reel0, renard_reel0, STEP, DAYS
    )
    mse_lapin = np.mean((lapin_pred - lapin_reel) ** 2)
    mse_renard = np.mean((renard_pred - renard_reel) ** 2)
    
    return mse_lapin + mse_renard

# Recherche des meilleurs paramètres
param_values = [1 / 3, 2 / 3, 1, 4 / 3]
best_params = [None]
best_mse = float('inf')

for alpha, beta, delta, gamma in product(param_values, repeat=4):
    mse_current = calculate_mse(alpha, beta, delta, gamma, lapin_reel[0], renard_reel[0])
    if mse_current < best_mse:
        best_mse = mse_current
        best_params = [alpha, beta, delta, gamma]

print("Meilleurs paramètres :", best_params)
print("Erreur MSE minimale :", best_mse)

# Simulation avec les meilleurs paramètres
lapin_pred, renard_pred, time = simulate_model_euler(
    best_params[0], best_params[1], best_params[2], best_params[3],
    lapin_reel[0], renard_reel[0], 
    STEP, DAYS*100
)

# Tracer les résultats
plt.figure(figsize=(15, 6))
plt.plot(time, lapin_pred, label="Lapins (prédit)", linewidth=2)
plt.plot(time, renard_pred, label="Renards (prédit)", linewidth=2)
plt.xlabel("Temps")
plt.ylabel("Population")
plt.title("Modèle proie-prédateur (Lapins et Renards)")
plt.legend()
plt.grid()
plt.show()
