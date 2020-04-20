from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["scikit-learn>=0.22", "numpy", "pandas>=1", "pysam>=0.15", "mappy>=2.17","matplotlib>=3", "seaborn>=0.10", "scipy"]

setup(
    name="pbampliconphasing",
    version="0.0.1",
    author="John Harting",
    author_email="jharting@pacifivbiosciences.com",
    description="Tools for clustering PB amplicons",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/PacificBiosciences/pbampliconclustering/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7"
    ],
)
