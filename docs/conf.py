# docs/source/conf.py

import os
import sys

# Make the project and Python API importable for autodoc.
project_root = os.path.abspath("../..")
python_api_path = os.path.join(project_root, "python_api")

sys.path.insert(0, project_root)
sys.path.insert(0, python_api_path)

project = "Met-ID"
author = "Your Name"
copyright = "2026, Your Name"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Uncomment if autodoc can't import dependencies in the Read the Docs build.
# Replace these with modules your Python API actually imports.
# autodoc_mock_imports = ["requests", "pandas", "psutil", "toml"]

autodoc_member_order = "bysource"
autodoc_typehints = "description"

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]