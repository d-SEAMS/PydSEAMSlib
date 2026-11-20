# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import importlib.util
import os
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, os.path.abspath("../../src"))
# dseams imports without its extension module for autodoc.
os.environ.setdefault("DSEAMS_DOC_BUILD", "1")
autodoc_mock_imports = ["pydseams.yoda"]

# Napoleon: NumPy docstrings in src/pydseams/*.py. yoda is mocked.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pydseams"
# the release string is the wheel version in pyproject.toml, so the docs
# cannot drift from the package
with (Path(__file__).resolve().parents[2] / "pyproject.toml").open("rb") as _fh:
    release = tomllib.load(_fh)["project"]["version"]
copyright = "2023\u2013present, d-SEAMS developers"
author = "Ruhila S (https://orcid.org/0000-0002-5443-9356)"


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx_contributors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_sitemap",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.mermaid",
]
autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = []

# The suffix(es) of source filenames.
# Authored pages are ox-rst from orgmode/; api.md and history.md stay Myst.
source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "shibuya"
html_title = "pydseams Python API documentation"
html_meta = {
    "description": (
        "Python APIs for d-SEAMS molecular simulation trajectory analysis, "
        "including ice, clathrate, ring, cage, and order-parameter methods."
    ),
}
html_static_path = ["_static"]
# Machine-readable discovery files at the site root.
html_extra_path = ["../../llms.txt", "../robots.txt"]
html_favicon = "_static/dseams/favicon.ico"
html_css_files = [
    "dseams/dseams-tokens.css",
    "dseams/dseams-family.css",
    "dseams/dseams-shibuya.css",
]

html_context = {
    "source_type": "github",
    "source_user": "d-SEAMS",
    "source_repo": "PydSEAMSlib",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

# Mermaid: vendored library; diagrams authorable via ``.. mermaid::`` (from Org RST export).
mermaid_version = "11.4.0"
mermaid_use_local = "dseams/mermaid.mjs"
mermaid_init_config = {
    "startOnLoad": False,
    "themeVariables": {"fontFamily": "ui-sans-serif, system-ui, sans-serif"},
}

html_theme_options = {
    "github_url": "https://github.com/d-SEAMS/PydSEAMSlib",
    "accent_color": "blue",
    "show_ai_links": False,
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    "toctree_collapse": True,
    "toctree_maxdepth": 3,
    "toctree_titles_only": True,
}

html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/edit-this-page.html",
    ],
}

html_baseurl = "https://python.dseams.info/"
sitemap_url_scheme = "{link}"

# --- Plugin options

myst_enable_extensions = [
    "deflist",
    "fieldlist",
    "substitution",
]

_names_path = Path(__file__).resolve().parents[2] / "scripts" / "yoda_public_names.py"
_names_spec = importlib.util.spec_from_file_location("yoda_public_names", _names_path)
if _names_spec is None or _names_spec.loader is None:
    raise RuntimeError(f"cannot load {_names_path}")
_names_mod = importlib.util.module_from_spec(_names_spec)
sys.modules[_names_spec.name] = _names_mod
_names_spec.loader.exec_module(_names_mod)
_yoda_names = _names_mod.public_names()
myst_substitutions = {
    "yoda_names": str(len(_yoda_names.exposed)),
    "yoda_functions": str(len(_yoda_names.exposed_functions)),
    "yoda_types": str(len(_yoda_names.exposed_types)),
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "ase": ("https://docs.ase-lib.org/", None),
    "dseams": ("https://docs.dseams.info", None),
    "luadseams": ("https://lua.dseams.info/", None),
}

bibtex_bibfiles = ["bibtex/pyseamsDocs.bib"]


def _org_file_links_to_html(app, doctree, docname):
    """Resolve the ``page.rst`` and ``page.md`` links that ox-rst writes for
    Org file links to the built pages; Sphinx treats them as external URIs."""
    from docutils import nodes

    srcdir = Path(app.srcdir).resolve()
    for ref in doctree.findall(nodes.reference):
        uri = ref.get("refuri", "")
        path, _, anchor = uri.partition("#")
        if "://" in path or not path.endswith((".rst", ".md")):
            continue
        target = (srcdir / os.path.dirname(docname) / path).resolve()
        if not target.is_file() or srcdir not in target.parents:
            continue
        target_doc = target.relative_to(srcdir).with_suffix("").as_posix()
        if target_doc not in app.env.all_docs:
            continue
        new = app.builder.get_relative_uri(docname, target_doc)
        ref["refuri"] = f"{new}#{anchor}" if anchor else new


def setup(app):
    app.connect("doctree-resolved", _org_file_links_to_html)
