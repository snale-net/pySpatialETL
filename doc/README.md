# Sphinx Documentation - SpatialETL

Technical guide for generating SpatialETL documentation using Sphinx.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Quick documentation generation](#quick-documentation-generation)
- [Documentation project structure](#documentation-project-structure)
- [Sphinx Configuration (`conf.py`)](#sphinx-configuration-confpy)
- [AutoAPI Templates Customization](#autoapi-templates-customization)
- [Project metadata](#project-metadata)
- [Theme and customization](#theme-and-customization)
- [Homepage structure (`index.rst`)](#homepage-structure-indexrst)
- [Useful commands](#useful-commands)
- [Resources](#resources)

## Prerequisites

- **UV installed** (see [UV documentation](https://docs.astral.sh/uv/))
- SpatialETL project cloned with UV environment configured

> **Note**: This guide uses UV commands (`uv run make html`). If you're not using UV,
> replace `uv run make` with standard `make` commands.

## Quick documentation generation

### 1️⃣ From project root
```bash
cd pySpatialETL
```

### 2️⃣ Navigate to `doc` folder
```bash
cd doc
```

### 3️⃣ Install requirements
```bash
uv run pip install -r requirements.txt
```

### 4️⃣ Generate HTML documentation
```bash
uv run make html
```

### 5️⃣ View the result

Open `doc/build/html/index.html` in your browser.

### 🔄 Complete rebuild

To force a full reconstruction:
```bash
uv run make clean && uv run make html
```

---

## Documentation project structure
```
doc/
├── source/             # Documentation source files
│   ├── conf.py         # Sphinx configuration (extensions, theme, paths)
│   ├── index.rst       # Documentation homepage
│   ├── api-reference/  # Auto api generate docs
│   ├── docs/           # User and developer documentation
│   ├── _static/        # Custom CSS files, images
│   └── _templates/     # Custom templates (autoapi, etc.)
├── build/              # Generated documentation (HTML, PDF, etc.)
│   └── html/           # HTML version of documentation
├── make.bat            # Build commands (Windows)
├── Makefile            # Build commands (Linux/macOS)
├── README.md           # This file
└── requirements.txt    # Sphinx dependencies (sphinx, autoapi, rtd-theme)
```

**Important**: The `build/` directory should **never** be committed to Git (see `.gitignore`).

---

## Sphinx Configuration (`conf.py`)

### Enabled extensions
```python
extensions = [
    'autoapi.extension',       # Automatically generate API documentation
    'sphinx.ext.autodoc',      # Extract docstrings
    'sphinx.ext.napoleon',     # Support Google/NumPy docstrings
    'sphinx.ext.viewcode',     # Links to source code
    'sphinx_rtd_theme',        # Read The Docs theme
    'sphinx_design',           # Design elements (buttons, grids)
]
```

### AutoAPI Configuration

AutoAPI automatically generates API documentation from source code.

**Standard configuration**:
```python
autoapi_generate_api_docs = True
autoapi_dirs = [
    # Core package
    "../../spatialetl-core/spatialetl",

    # Common providers
    "../../providers/common/gdal/spatialetl",
    "../../providers/common/grib/spatialetl",
    "../../providers/common/mpi/spatialetl",
    "../../providers/common/netcdf/spatialetl",

    # Data source providers
    "../../providers/ecmwf/spatialetl",
    "../../providers/gmt/spatialetl",
    "../../providers/hycom/spatialetl",
    "../../providers/mercator/spatialetl",
    "../../providers/meteofrance/spatialetl",
    "../../providers/swan/spatialetl",
    "../../providers/symphonie/spatialetl",
    "../../providers/telemac/spatialetl",
    "../../providers/ww3/spatialetl",
]
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

### File exclusion
```python
exclude_patterns = [
    '_build',
    'build',
    '**/tests/*'
]
```

These patterns prevent Sphinx from processing generated files or tests.

---

## AutoAPI Templates Customization

The project uses custom AutoAPI templates located in `source/_templates/autoapi/`.

### Template structure
```
_templates/autoapi/
├── index.rst           # Main API reference index page
├── macros.rst          # Reusable Jinja2 macros
└── python/             # Python-specific templates
    ├── attribute.rst   # Attribute documentation
    ├── class.rst       # Class documentation
    ├── data.rst        # Data/constants documentation
    ├── exception.rst   # Exception documentation
    ├── function.rst    # Function documentation
    ├── method.rst      # Method documentation
    ├── module.rst      # Module documentation
    ├── package.rst     # Package documentation
    └── property.rst    # Property documentation
```

### How it works

AutoAPI uses these templates to generate documentation in `source/api_reference/`.
The templates use Jinja2 syntax and can access object metadata (names, docstrings, signatures, etc.).

### Modifying templates

To customize the API documentation appearance:
1. Edit the relevant `.rst` template in `source/_templates/autoapi/python/`
2. Rebuild documentation: `uv run make clean && uv run make html`
3. Check the result in `build/html/api_reference/`

---

## Project metadata

In `conf.py`:
```python
from datetime import datetime

project = 'SpatialETL documentation'
copyright = f'{datetime.now().year}, SNALE'
author = 'SNALE'
release = '0.1'
```

---

## Theme and customization

### Read The Docs Theme
```python
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
```

The `source/_static/custom.css` file allows adding custom styles.

---

## Homepage structure (`index.rst`)

The `source/index.rst` file defines the navigation structure.

---

## Useful commands

| Command | Description |
|---------|-------------|
| `uv run make html` | Generate HTML documentation |
| `uv run make clean` | Remove generated files |
| `uv run make clean && uv run make html` | Clean and regenerate |
| `uv run make linkcheck` | Check external links |

## Resources

- [Official Sphinx Documentation](https://www.sphinx-doc.org/)
- [Sphinx AutoAPI](https://sphinx-autoapi.readthedocs.io/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [Read The Docs Theme](https://sphinx-rtd-theme.readthedocs.io/)
- [UV Documentation](https://docs.astral.sh/uv/)