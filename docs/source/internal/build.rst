Build Process
*************

Github Actions
==============
ARC uses Github actions to perform a variety of useful build steps

**Testing (CI/CD)**

ARC has a series of tests that ensure that the project is not FUBAR.

These consist of some basic unit tests
- Ensuring output is correct
- Ensuring invalid characters throw a useful error
- Ensuring that empty files throw an error

These are useful when you are uploading changes so there is a Github Actions pipeline set up

The following code for automated testing is found in the repo as ``.github/workflows/ARC-test.yaml``

.. code-block:: yaml

    on:
    push:
        branches: [ master ]
    pull_request:
        branches: [ master ]

This tells Github actions to run this workflow when something is pushed to master or when a pull request is made for master

.. code-block:: yaml

    jobs:
    build:
        # Specify we want to run on the latest ubuntu image
        runs-on: ubuntu-latest

        steps:
        - uses: actions/checkout@v2
        - name: Set up Python 3.9
        uses: actions/setup-python@v2
        with:
            python-version: 3.9
        # Run a series of commands to install dependencies (Ubuntu apt-get install command)
        - name: Install dependencies
        run: |
            sudo apt-get install hmmer ncbi-blast+ git
            python -m pip install --upgrade pip
            pip install flake8 pytest
            if [ -f requirements.txt ]; then pip install -r requirements.txt; fi
        # Linting just ensures that there are no syntax or formatting errors
        - name: Lint with flake8
        run: |
            # stop the build if there are Python syntax errors or undefined names
            flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
            # exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
            flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
        # Actually run the tests with Pytest
        - name: Test with pytest
        run: |
            pytest

Pytest can be run from the main directory, it scans inside the test folder and will autofind tests/arc_test.py

If any of the 4 tests fail the build will fail and this is represented as an X next to your commit on the top window in Github

You can also view the status and a series of logs for Github Actions under the "Actions" tab on the repository (https://github.com/IEDB/ARC/actions/runs/951740407) as an example

Development Tips
================

Please ensure that when you are making a new release on Github you properly iterate the versioning in setup.py in the top level of the directory

.. code-block:: python
    
    import setuptools

    with open("README.md", "r") as fh:
        long_description = fh.read()

    setuptools.setup(
        name="bio-arc",
        # Right here
        version="0.1.1",

When you are installing the package, install using ``pip install -e .`` from the top directory. This lets you edit the package without reinstalling each time.

Similarly, use a virtualenv when you are developing to keep things clean.