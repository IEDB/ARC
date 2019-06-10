import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bio-arc",
    version="0.0.1",
    author="Austin Crinklaw",
    author_email="acrinklaw@lji.org",
    description="Antigen Receptor Classifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iedb/arc",
    packages=setuptools.find_packages(),
    package_data={'ARC': ['data/HMMs/*', 'data/MHC_HMMs/*', 'data/MRO_Gdomain.csv', 'data/blastdb/*', 'build_pipeline/*']},
    install_requires=[
        'pandas',
        'biopython',],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: Unix",
    ],
)
