# Conf.py Configuration file for the Sphinx documentation builder.

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../.."))
import spatialetl

project = 'pySpatialETL documentation'
copyright = f'{datetime.now().year}, SNALE'
author = 'SNALE'
release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    'autoapi.extension', # Génère automatiquement une API doc depuis le code source
    'sphinx.ext.autodoc',  # Génération automatique de la doc depuis les docstrings
    'sphinx.ext.napoleon',  # Support des formats Google et NumPy dans les docstrings
    'sphinx.ext.viewcode', # Ajoute des liens vers le code source
    'sphinx_rtd_theme',
    'sphinx.ext.autosummary',
    'sphinx_design',
]

templates_path = ['_templates']

exclude_patterns = []

autosummary_generate = False

# -- AutoAPI configuration ---------------------------------------------------
autoapi_generate_api_docs = True
autoapi_dirs = ["../../spatialetl"]
autoapi_root = "api_reference"
autoapi_keep_files = True            # Conserve les fichiers .rst générés
autoapi_add_toctree_entry = False    # Ne pas injecter automatiquement dans la toctree (on le gère à la main)

autoapi_options = [
    "members", # "members": inclut les membres (attributs et méthodes) des classes.
    "undoc-members", # "undoc-members": inclut les membres non documentés.
    "private-members", # "private-members": inclut les membres privés (commençant par un underscore).
    "show-inheritance", # "show-inheritance": affiche l'héritage des classes.
    "show-module-summary", # "show-module-summary": affiche un résumé du module.
    "special-members", # "special-members": inclut les méthodes spéciales (par exemple, __init__, __str__).
]

autoapi_ignore = [
    # "**/Coverage.py",
]

# -- Napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme' # Thème ReadTheDocs
html_static_path = ['_static'] # Dossier pour les CSS/images statiques
html_css_files = [
    'custom.css',
]