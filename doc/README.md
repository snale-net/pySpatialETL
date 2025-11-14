# Sphinx Documentation - pySpatialETL

Technical guide for generating pySpatialETL documentation using Sphinx.

## Prerequisites

- **Python 3.9+**
- **UV installed** (see [UV documentation](https://docs.astral.sh/uv/))
- pySpatialETL project cloned with UV environment configured

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

### 3️⃣ Generate HTML documentation
```bash
uv run make html
```

### 4️⃣ View the result

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
├── source/              # Documentation source files
│   ├── conf.py         # Sphinx configuration (extensions, theme, paths)
│   ├── index.rst       # Documentation homepage
│   ├── docs/           # User and developer documentation
│   ├── _static/        # Custom CSS files, images
│   └── _templates/     # Custom templates (autoapi, etc.)
├── build/              # Generated documentation (HTML, PDF, etc.)
│   └── html/           # HTML version of documentation
├── Makefile            # Build commands (Linux/macOS)
├── make.bat            # Build commands (Windows)
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
autoapi_dirs = ["../../spatialetl-core/spatialetl"]  # UV workspace architecture
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

**Reference**: [AutoAPI Template Documentation](https://sphinx-autoapi.readthedocs.io/en/latest/reference/templates.html)

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

### Important directives

- **toctree**: Defines the table of contents
- **:maxdepth:**: Display depth of tree structure
- **:caption:**: Section title

---

## Useful commands

| Command | Description |
|---------|-------------|
| `uv run make html` | Generate HTML documentation |
| `uv run make clean` | Remove generated files |
| `uv run make clean && uv run make html` | Clean and regenerate |
| `uv run make linkcheck` | Check external links |

---

## Troubleshooting

### "ModuleNotFoundError" error

**Symptom**: Sphinx cannot find Python modules to document.

**Solution**: Verify UV environment is synchronized:
```bash
cd pySpatialETL
uv sync
```

### Warnings "document not in toctree"

**Symptom**: `.rst` files created but not referenced in a `.. toctree::`.

**Solutions**:
1. Add the file to a `toctree` block in `index.rst` or a parent file
2. Add the file pattern to `exclude_patterns` in `conf.py` if it shouldn't be included

---

## Best practices

### ✅ Do

- Always use `uv run make html` to ensure correct environment usage
- Test documentation locally before committing sources
- Follow reStructuredText conventions for `.rst` files
- Document functions/classes with Google or NumPy docstrings

### ❌ Don't

- **Never** commit the `build/` directory
- Don't manually modify files generated by autoapi in `api_reference/`
- Avoid absolute paths in `.rst` files

---

## Resources

- [Official Sphinx Documentation](https://www.sphinx-doc.org/)
- [Sphinx AutoAPI](https://sphinx-autoapi.readthedocs.io/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [Read The Docs Theme](https://sphinx-rtd-theme.readthedocs.io/)
- [UV Documentation](https://docs.astral.sh/uv/)

---