import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bio-arc",
    version="0.1.0",
    author="Austin Crinklaw",
    author_email="acrinklaw@lji.org",
    description="Antigen Receptor Classifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iedb/arc",
    packages=setuptools.find_packages(),
    package_data={
        'ARC': [
            'data/*', 'data/blastdb', 'data/HMMs', 'data/IgNAR',
            'data/MHC_HMMs', 'tests/*'
        ]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
)
