import os
import sys
sys.path.insert(0, os.path.abspath('../../'))  # project root

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Custom Bioinformatics Tools'
copyright = '2026, Eric Moeller'
author = 'Eric Moeller'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',      # Core extension to extract docstrings
    'sphinx.ext.napoleon',     # Needed if using Google or NumPy styles
    'sphinx.ext.viewcode'      # Optional: Adds links to highlighted source code
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
