import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../.."))


project = 'pySpatialETL documentation'
copyright = f'{datetime.now().year}, SNALE'
author = 'SNALE'
release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    'autoapi.extension',      # Automatically generate API docs from source code
    'sphinx.ext.autodoc',     # Auto-generate documentation from docstrings
    'sphinx.ext.napoleon',    # Support for Google and NumPy docstring formats
    'sphinx.ext.viewcode',    # Add links to highlighted source code
    'sphinx_rtd_theme',       # ReadTheDocs theme for HTML output
    'sphinx.ext.autosummary', # Generate summary tables for modules/classes/functions
    'sphinx_design',          # Enhanced design elements (buttons, grids, etc.)
]

autodoc_mock_imports = [
    "numpy",
    "scipy",
    "pandas",
    "array_split",
#setuptools, cython
]
templates_path = ['_templates']

exclude_patterns = []

autosummary_generate = True

# -- AutoAPI configuration ---------------------------------------------------
autoapi_generate_api_docs = True
autoapi_dirs = ["../../spatialetl"]
autoapi_root = "api_reference"
autoapi_keep_files = True            # Keep generated .rst files
autoapi_add_toctree_entry = False    # Do not auto-inject into the toctree (manual control)

autoapi_options = [
    "members",              # Include class members (attributes and methods)
    "undoc-members",        # Include undocumented members
    "private-members",      # Include private members (starting with _)
    "show-inheritance",     # Show class inheritance
    "show-module-summary",  # Show a summary at the top of each module
    "special-members",      # Include special methods (e.g., __init__, __str__)
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
html_theme = 'sphinx_rtd_theme'            # Use the ReadTheDocs HTML theme
html_static_path = ['_static']             # Directory for static CSS/images
html_css_files = [
    'custom.css',
]