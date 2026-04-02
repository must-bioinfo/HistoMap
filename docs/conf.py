# docs/conf.py
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# Project information
project = 'HistoMapTx'
copyright = '2023-2025, HistoMapTx Contributors'
author = 'HistoMapTx Contributors'
version = '0.1.0'
release = '0.1.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',       # API docs from docstrings
    'sphinx.ext.viewcode',      # View source code
    'sphinx.ext.napoleon',      # Support NumPy and Google docstrings
    'sphinx.ext.intersphinx',   # Link to other projects' docs
    'sphinx_autodoc_typehints', # Type hints in docs
    'nbsphinx',                 # Jupyter notebook support
    'sphinx_copybutton',        # Add copy buttons to code blocks
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = False
napoleon_use_ivar = True

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'special-members': '__init__',
}

# Theme settings
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']
html_logo = '../figures/logo.png'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'navigation_depth': 4,
}

# Notebook settings
nbsphinx_execute = 'never'
nbsphinx_allow_errors = True
nbsphinx_timeout = 600

# Cross-reference other documentation
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
}

# Ignore certain files/patterns
exclude_patterns = ['_build', '**.ipynb_checkpoints']