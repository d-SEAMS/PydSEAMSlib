# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

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
copyright = "2023--present, d-SEAMS developers"
author = "Ruhila S"


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
# llms.txt at the site root, from the repository root.
html_extra_path = ["../../llms.txt"]
html_favicon = "_static/logo/pydseams-icon.ico"
html_css_files = ["custom.css"]

html_context = {
    "source_type": "github",
    "source_user": "d-SEAMS",
    "source_repo": "PydSEAMSlib",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

# Mermaid: use default CDN; diagrams authorable via ``.. mermaid::`` (from Org RST export).
mermaid_version = "11.4.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true, theme:'neutral'});"

html_theme_options = {
    "github_url": "https://github.com/d-SEAMS/PydSEAMSlib",
    "accent_color": "teal",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    "toctree_collapse": True,
    "toctree_maxdepth": 3,
    "toctree_titles_only": True,
    "light_logo": "_static/logo/pydseams-logo-light.png",
    "dark_logo": "_static/logo/pydseams-logo-dark.png",
    "nav_links": [
        {
            "title": "Ecosystem",
            "children": [
                {
                    "title": "d-SEAMS engine",
                    "url": "https://docs.dseams.info",
                    "summary": "libyodaLib and the seams CLI",
                    "external": True,
                },
                {
                    "title": "pydseams",
                    "url": "https://python.dseams.info/",
                    "summary": "Python Frame API on yoda",
                    "external": True,
                },
                {
                    "title": "dseams (Lua)",
                    "url": "https://lua.dseams.info/",
                    "summary": 'require("dseams") and Fennel',
                    "external": True,
                },
                {
                    "title": "linkcell",
                    "url": "https://github.com/d-SEAMS/linkcell",
                    "summary": "Periodic linked-cell k-nearest neighbours",
                    "external": True,
                },
            ],
        },
    ],
}

html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
    ],
}

html_baseurl = "https://python.dseams.info/"
sitemap_url_scheme = "{link}"

# --- Plugin options

myst_enable_extensions = [
    "deflist",
    "fieldlist",
]

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
