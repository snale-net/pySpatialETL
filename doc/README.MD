# Sphinx Installation and Documentation Deployment V1

## Requirements

- **Python 3.9**

## 1️⃣ Create and Activate a Virtual Environment

```sh
python -m venv .venv
source .venv/bin/activate
```

## 2️⃣ Navigate to the `doc` Folder 📂

```sh
cd doc
```

## 3️⃣ Install Sphinx and Required Extensions

```sh
pip install -r requirements.txt
```

The following extensions are used:

- `sphinx==7.4`
- `sphinx-autoapi~=3.6` → Automatically generates documentation from source code.
- `sphinx-rtd-theme~=3.0` → Applies the **Read the Docs** theme, widely used for Python documentation.

## 4️⃣ Generate Documentation

### 📌 **Build the Documentation in HTML**

Inside the `doc` folder, run:

```sh
make html
```

## 🎯 Once the generation is complete, open `doc/build/html/index.html` in a browser to view the documentation.

To regenerate it, use:

```sh
make clean && make html
```

---

## 1️⃣PS: Initial Setup

### 🛠️ Configuration Command (from `doc` folder)

```sh
sphinx-quickstart
```

Configuration Steps:

- **Separate source and build directories?** → **Y** (Yes)
- **Project name:** → `pySpatialETL`
- **Author name(s):** → `SNALE`
- **Project version:** → `0.1`
- **Project language:** → (leave empty, defaults to `en`)

After this setup, Sphinx generates the following structure:

```
📂 doc/                # Root folder for documentation
├── 📂 source/         # Contains documentation source files
│   ├── 📄 conf.py        # Main configuration file (project settings, extensions, etc.)
│   ├── 📄 index.rst      # Main index file (homepage of documentation)
│   ├── 📂 _static/       # Folder for static files (CSS, images, etc.)
│   └── 📂 _templates/    # Folder for custom templates
├── 📂 build/          # Stores generated documentation (HTML, LaTeX, etc.)
├── 📄 Makefile        # Build file for Linux/macOS
└── 📄 make.bat        # Build file for Windows
```

## 2️⃣ Define Project Information

Inside `conf.py`, define project metadata:

```python
from datetime import datetime
project = 'pySpatialETL'
copyright = f'{datetime.now().year}, SNALE'
author = 'SNALE'
release = '0.1'
```

---

## 3️⃣ Configure Project Path

```python
import os
import sys
sys.path.insert(0, os.path.abspath("../.."))
```

This adds the project's root directory to `sys.path`, allowing Sphinx to locate the source code for documentation.

---

## 4️⃣ Enabled Extensions

The project uses multiple extensions to enhance documentation:

```python
extensions = [
    'autoapi.extension',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
    'sphinx.ext.autosummary',
]
```

- **autoapi.extension**: Automatically generates API documentation.
- **sphinx.ext.autodoc**: Extracts docstrings from Python code.
- **sphinx.ext.napoleon**: Supports Google and NumPy docstring formats.
- **sphinx.ext.viewcode**: Adds links to the source code within the documentation.
- **sphinx_rtd_theme**: Applies the Read The Docs theme.
- **sphinx.ext.autosummary**: Generates automatic summaries of modules and classes.

---

## 5️⃣ AutoAPI Configuration

AutoAPI is used to automatically generate module documentation:

```python
autoapi_generate_api_docs = True
autoapi_dirs = ["../../spatialetl"]
autoapi_root = "api_reference"
autoapi_keep_files = True
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
]
```

- **autoapi_generate_api_docs**: Enables automatic documentation generation.
- **autoapi_dirs**: Defines the directory containing source code to be documented.
- **autoapi_root**: Directory where the API documentation will be generated.
- **autoapi_keep_files**: Keeps generated files for review.
- **autoapi_options**: Customizes the documentation output.

---

## 6️⃣ Theme and Static File Configuration

```python
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
```

- **html_theme**: Applies the Read The Docs theme.
- **html_static_path**: Directory for static files.

---

## 7️⃣ Index File Structure (`index.rst`)

The `index.rst` file serves as the homepage and includes:

```
Welcome to pySpatialETL documentation!
======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

- **toctree**: Defines the table of contents.
- **genindex**: Adds a general index.
- **modindex**: Adds a module index.
- **search**: Enables a search function.

---

