import sys
import os

# Absolute path to the root of your project, relative to this conf.py
conf_dir = os.path.dirname(os.path.abspath(__file__))  # docs/source
project_root = os.path.abspath(os.path.join(conf_dir, "..", ".."))
sys.path.insert(0, project_root)

# Mock heavy scientific packages for autodoc
autodoc_mock_imports = [
    "obsplus",
    "obspy",
    "datasets",
    "huggingface_hub",
    "pyarrow",
    "matplotlib",
    "seaborn",
    "scipy",
    "pyarrow",
    "pandas",
    "gamma",
    "EQTransformer",
]

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SeisMonitor'
copyright = '2025, Emmanuel Castillo'
author = 'Emmanuel Castillo'
release = '0.0.58'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#extensions = []
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "myst_parser",
    "sphinx_autodoc_typehints",
    "sphinx_design"
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


rst_prolog = """
.. role:: red
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

html_theme_options = {
    'logo_only': True,
    'display_version': False
}
html_logo = "_static/logo/seismonitor.PNG"

def setup(app):
    app.add_css_file('css/custom.css')
