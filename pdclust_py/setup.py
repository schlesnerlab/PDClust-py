from setuptools import setup, find_packages

setup(
    name='pdclust_py',
    version='0.0.1',
    description='Python reimplementation of the clustering approach used in "High-Resolution Single-Cell DNA Methylation Measurements Reveal Epigenetically Distinct Hematopoietic Stem Cell Subpopulations"',
    author='Theodor F. Pöschl',
    url='https://github.com/TeaBag-T/pdclust_py',
    packages=find_packages(),
    keywords='pdclust clustering single-cell pairwise dissimilarity',
    install_requires=[
        'numpy',
        'pandas',
        'seaborn',
        'PyComplexHeatmap',
        'dask',
    ],
)
