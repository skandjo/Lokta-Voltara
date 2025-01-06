# README

## Description
Ce projet implémente un modèle proie-prédateur basé sur les équations de Lotka-Volterra. Le fichier principal, **`main.py`**, permet de simuler l'évolution des populations de lapins et de renards sur une période donnée.

## Instructions d'exécution
1. Placez tous les fichiers nécessaires dans le même répertoire (notamment **`main.py`** et **`populations_lapins_renards.csv`**).
2. Exécutez le script avec la commande :
   ```bash
   python main.py
   ```
3. Le script génère des graphiques et affiche les meilleurs paramètres ainsi que l'erreur MSE dans la console.

## Remarques
- Le choix du nombre de jours (`days`) et du pas (`step`) est crucial pour obtenir des résultats précis.
- Le MSE et le Grid Search fonctionnent correctement, mais les valeurs optimales de `days` et `step` doivent être expérimentées.

