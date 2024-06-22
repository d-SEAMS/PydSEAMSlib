# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pyseams'
copyright = '2024, RuhiRG'
author = 'RuhiRG'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#referenced my favourite [1]
extensions = [
    "myst_parser",
    "autodoc2",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx_contributors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.spelling",
]

autodoc2_render_plugin = "myst"
autodoc2_packages = [
    f"../../src/{project}",
]

myst_enable_extensions = [
    "deflist",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']

#references
# [1] https://github.com/HaoZeke/openblas_buildsys_snips/blob/main/docs/source/conf.py