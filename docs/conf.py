#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
author = "PrOMMiS Contributors"
autosummary_generate = True
comments_config = {"hypothesis": False, "utterances": False}
copyright = "2026"
exclude_patterns = ["**.ipynb_checkpoints", ".DS_Store", "Thumbs.db", "_build"]
extensions = [
    "sphinx_togglebutton",
    "sphinx_copybutton",
    "myst_nb",
    "jupyter_book",
    "sphinx_thebe",
    "sphinx_comments",
    "sphinx_external_toc",
    "sphinx.ext.intersphinx",
    "sphinx_design",
    "sphinx_book_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "prommis._sphinx_extensions.chem",
    "sphinx_jupyterbook_latex",
    "sphinx_multitoc_numbering",
]
external_toc_exclude_missing = False
external_toc_path = "_toc.yml"
html_baseurl = ""
html_favicon = ""
html_logo = "images/PROMMIS Logo_Stacked Color_GreyText_TransparentBackground.png"
html_sourcelink_suffix = ""
html_theme = "sphinx_book_theme"
html_theme_options = {
    "search_bar_text": "Search this book...",
    "launch_buttons": {
        "notebook_interface": "classic",
        "binderhub_url": "",
        "jupyterhub_url": "",
        "thebe": False,
        "colab_url": "",
        "deepnote_url": "",
    },
    "path_to_docs": "",
    "repository_url": "https://github.com/executablebooks/jupyter-book",
    "repository_branch": "master",
    "extra_footer": "",
    "home_page_in_toc": True,
    "announcement": "",
    "analytics": {
        "google_analytics_id": "",
        "plausible_analytics_domain": "",
        "plausible_analytics_url": "https://plausible.io/js/script.js",
    },
    "use_repository_button": False,
    "use_edit_page_button": False,
    "use_issues_button": False,
}
html_title = "PrOMMiS"
latex_engine = "pdflatex"
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_admonition",
    "linkify",
    "substitution",
    "tasklist",
    "dollarmath",
]
myst_url_schemes = ["mailto", "http", "https"]
nb_execution_allow_errors = False
nb_execution_cache_path = ""
nb_execution_excludepatterns = []
nb_execution_in_temp = False
nb_execution_mode = "off"
nb_execution_timeout = 30
nb_output_stderr = "show"
numfig = True
pygments_style = "sphinx"
suppress_warnings = ["etoc.toctree"]
templates_path = ["_templates"]
use_jupyterbook_latex = True
use_multitoc_numbering = True
