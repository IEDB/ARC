API Documenation
****************

Generally speaking other code will create a Classifier object that will be able to accept a protein sequence and return a classification result

Example
=======

Here is a segment of code used in the PDB Classifier

.. code-block:: python

    from ARC.classifier import SeqClassifier

    classifier = SeqClassifier()
    # entry[0] refers to ID, entry[1] refers to sequence
    res = classifier.gen_classify(entry[1], entry[0])

This will return a tuple ``(receptor, chain_type, calc_mhc_allele)``

Further documentation for the classifier.py and mhc_G_domain.py files can be found below

Classifier
==========
.. automodule:: ARC.classifier
    :members:

MHC G Domain
============
.. automodule:: ARC.mhc_G_domain
    :members:
