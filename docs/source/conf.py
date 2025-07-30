import os
import sys
# sys.path.insert(0, os.path.abspath('.')) # This points to the current dir (docs/source), not your project root
sys.path.insert(0, os.path.abspath('../../src')) # Adjust this path to point to your project root
# If your source code is directly in 'my_project_root/my_package_name', you might need:
# sys.path.insert(0, os.path.abspath('../..'))
# If your source code is in 'my_project_root/src/my_package_name', you might need:
# sys.path.insert(0, os.path.abspath('../../src'))
# 
# 
# # Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DNRlib'
copyright = '2025, CITCEA'
author = 'CITCEA'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon', # If using Google/NumPy style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx_autodoc_typehints', # For better type hint rendering
]

# Optional: Configure intersphinx for linking to other projects' docs
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'requests': ('https://requests.readthedocs.io/en/latest/', None),
    # Add more mappings as needed
}

# Optional: Configure Napoleon (for Google/NumPy docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = False # Set to True if using NumPy style exclusively
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_args = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Optional: Configure Type Hints
autodoc_typehints = "signature" # "signature" to put type hints in the function signature, "description" to put them in the description section
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False, # Set to True to include private members (e.g., _private_method)
    'special-members': False, # Set to True to include special members (e.g., __init__)
    'inherited-members': False,
    'show-inheritance': True,
}

# Set the theme (if you installed sphinx-rtd-theme)
html_theme = 'sphinx_rtd_theme'

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_static_path = ['_static']
