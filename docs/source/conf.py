# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path
# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../'))
HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / 'extensions')]

on_rtd = os.environ.get('READTHEDOCS') == 'True'


# -- Project information -----------------------------------------------------

project = 'gefpy'
copyright = '2022, Huang Zhibo'
author = 'Huang Zhibo'

# The full version, including alpha/beta/rc tags
release = '0.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
#    'sphinx.ext.intersphinx',
#    'sphinx.ext.doctest',
#    'sphinx.ext.coverage',
#    'sphinx.ext.mathjax',
#    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
#    'nbsphinx',
#    "sphinx.ext.viewcode",
#    "sphinx_autodoc_typehints",
#    *[p.stem for p in (HERE / 'extensions').glob('*.py')],
]
#
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_default_flags = ['members']
todo_include_todos = False
api_dir = HERE / 'api'  # function_images
print(api_dir)
# The master toctree document.
language = None

#html_theme = 'sphinx_rtd_theme'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'default'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
