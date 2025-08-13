import os
import sys
import toml
sys.path.insert(0, os.path.abspath('..'))

project = 'pybhpt'
copyright = '2025, znasipak'
author = 'znasipak'

# Read version from pyproject.toml
pyproject_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml'))
with open(pyproject_path, 'r') as f:
    pyproject = toml.load(f)
release = pyproject['project']['version']

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'myst_nb',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',]
napoleon_use_ivar = True
myst_enable_extensions = ["dollarmath", "amsmath"]
nb_execution_mode = "off"
autosummary_generate = True
myst_dmath_double_inline = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__call__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_theme_options = {
    "repository_url": "https://github.com/znasipak/pybhpt",
    "use_repository_button": True,
    "navbar_end": ["navbar-icon-links"],
}
html_title = "pybhpt"
html_context = {
    "default_mode": "light"
}
