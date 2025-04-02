# 📌 Installing Sphinx for pySpatialETL Documentation  

## Prerequisites  

1. **Make sure you have Python 3.12 installed**:  
[🔗 Download Python 3.12](https://www.python.org/downloads/release/python-3120/)  

2. **Create and activate a virtual environment**  

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On macOS/Linux
   .venv\Scripts\activate     # On Windows
   ```  

3. **Install Sphinx**:  
   ```bash
   pip install sphinx==8.2.3
   ```  
   or  

   ```bash
   pip install -r doc/requirements.txt
   ```  

---

## Preconfigured Documentation Environment  
⚠️ **You do not need to perform the steps below, as they have already been completed.**  

1. **Create and activate a virtual environment**  

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On macOS/Linux
   .venv\Scripts\activate     # On Windows
   ```  

2. **Install Sphinx**  

   ```bash
   pip install sphinx==8.2.3
   ```  
   or  

   ```bash
   pip install -r doc/requirements.txt
   ```  

3. **Navigate to the `doc` folder 📂**  

   ```bash
   cd doc
   ```

4. **Initialize Sphinx**  

   ```bash
   sphinx-quickstart
   ```

   **Recommended answers**:  
   - **Separate source and build directories?** → **Y** (Yes)  
   - **Project name:** → `pySpatialETL`  
   - **Author(s) name(s):** → `SNALE`  
   - **Project version:** → `0.1`  
   - **Project language:** → *(leave empty, default is `en`)*  

---

## 3️⃣ Documentation Structure  

Once Sphinx has been initialized following the instructions above, the structure will look like this:  

```
📂 doc/                # Root folder for the documentation
├── 📂 source/         # Contains documentation source files
│   ├── 📄 conf.py        # Main configuration file (project settings, extensions, etc.)
│   ├── 📄 index.rst      # Main index file (documentation homepage)
│   ├── 📂 _static/       # Folder for static files (CSS, images, etc.)
│   └── 📂 _templates/    # Folder for custom templates
├── 📂 build/          # Stores the generated documentation (HTML, LaTeX, etc.)
├── 📄 Makefile        # File to generate documentation on Linux/macOS
└── 📄 make.bat        # File to generate documentation on Windows
```

## 4️⃣ Dynamic Year Replacement conf.py
- **`{{ current_year }}`**: This placeholder is meant to be dynamically replaced with the current year when you generate or update your documentation.

```python
from datetime import datetime

current_year = datetime.now().year
readme_content = f'{current_year}, SNALE'
```