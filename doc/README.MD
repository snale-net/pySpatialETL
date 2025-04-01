# 📌 Installation de Sphinx pour la documentation de pySpatialETL  

##  Prérequis  

1. **Assurez-vous d'avoir installé Python 3.12** :  
[🔗 Télécharger Python 3.12](https://www.python.org/downloads/release/python-3120/)  

2. **Créer et activer un environnement virtuel**  

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Sur macOS/Linux
   .venv\Scripts\activate    # Sur Windows
   ```  

3. **Installer Sphinx** :
   ```bash
   pip install sphinx==8.2.3
   ```  
    ou

   ```bash
   pip install -r doc/requirements.txt
   ```  

## Environnement de documentation préconfiguré
⚠️ **Vous n’avez pas besoin de faire les étapes décrites ci-dessous, elles ont déjà été effectuées.**  

1. **Créer et activer un environnement virtuel**  

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Sur macOS/Linux
   .venv\Scripts\activate    # Sur Windows
   ```  

2. **Installer Sphinx**  

   ```bash
   pip install sphinx==8.2.3
   ```  
    ou

   ```bash
   pip install -r doc/requirements.txt
   ```  

3. **Se déplacer dans le dossier `doc` 📂**  

   ```bash
   cd doc
   ```

4. **Initialiser Sphinx**  

   ```bash
   sphinx-quickstart
   ```

   **Réponses recommandées** :  
   - **Séparer les répertoires source et de sortie ?** → **Y** (Oui)  
   - **Nom du projet :** → `pySpatialETL documentation`  
   - **Nom(s) de(s) l’auteur(s) :** → `SNALE TEAM`  
   - **Version du projet :** → `0.1`  
   - **Langue du projet :** → *(laisser vide, par défaut `en`)*  

## 3️⃣ Structure de la documentation  

Une fois Sphinx initialisé en suivant les consignes ci-dessus, la structure ressemble à ceci :  

```
📂 doc/                # Dossier racine pour la documentation
├── 📂 source/         # Contient les fichiers source de la documentation
│   ├── 📄 conf.py        # Fichier de configuration principal (paramètres du projet, extensions, etc.)
│   ├── 📄 index.rst      # Fichier d'index principal (page d'accueil de la documentation)
│   ├── 📂 _static/       # Dossier pour les fichiers statiques (CSS, images, etc.)
│   └── 📂 _templates/    # Dossier pour les templates personnalisés
├── 📂 build/          # Stocke la documentation générée (HTML, LaTeX, etc.)
├── 📄 Makefile        # Fichier pour générer la documentation sur Linux/macOS
└── 📄 make.bat        # Fichier pour générer la documentation sur Windows
```
