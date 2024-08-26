# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "PydSEAMSlib"
release = "0.0.2"
copyright = "2024, d-SEAMS developers"
author = "Ruhila"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# referenced my favourite [1]
extensions = [
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx_contributors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
]
autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = []

# The suffix(es) of source filenames.
source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_title = "PydSEAMSlib"
html_static_path = ["_static"]

html_theme_options = {
    "repository_url": "https://github.com/d-SEAMS/PydSEAMSlib",
    "use_repository_button": True,
    "logo": {
        "image_light": "_static/logo/pydseamslib_logo_light.png",
        "image_dark": "_static/logo/pydseamslib_logo_dark.png",
    },
}
# --- Plugin options

myst_enable_extensions = [
    "deflist",
    "fieldlist",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

bibtex_bibfiles = ["bibtex/pyseamsDocs.bib"]

# references
# [1] https://github.com/HaoZeke/openblas_buildsys_snips/blob/main/docs/source/conf.py
