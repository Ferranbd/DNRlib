# setup.py
from setuptools import setup, find_packages

setup(
    name='DNRlib', # Replace with your actual project name
    version='0.1.0',          # Replace with your project's version
    packages=find_packages(where='.', exclude=('docs', 'tests')), # Adjust 'where' if your code is in 'src'
    # Add other metadata like:
    # install_requires=[
    #     'requests',
    #     'numpy',
    # ],
    description='Distribution Network Reconfiguration methods library',
    long_description=open('README.rst').read(),
    # long_description_content_type='text/markdown',
    # url='http://your-project-url.com',
    author='Ferran Bohigas i Daranas / CITCEA / UPC',
    author_email='ferran.bohigas@upc.edu',
)