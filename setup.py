from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()


setup(name = "allcopol",
    version = "0.1.1",
    description = "AllCoPol: Inferring allele co-ancestry in polyploids",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/AGOberprieler/allcopol",
    author = "Ulrich Lautenschlager",
    author_email = "ulrich.lautenschlager@ur.de",
    license = "MIT",
    packages = find_packages(),
    install_requires = [ 
        "argparse", "biopython", "configargparse", "numpy", "scipy"
    ],
    entry_points = {
        "console_scripts": [
            "allcopol=allcopol.allcopol:main",
            "align_clusters=allcopol.align_clusters:main",
            "create_indfile=allcopol.create_indfile:main",
            "relabel_trees=allcopol.relabel_trees:main",
        ],
    },
    zip_safe = False,
    python_requires = ">=3.5",
)

