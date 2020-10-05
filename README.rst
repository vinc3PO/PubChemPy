PubChemPy
=========

Fork of the glorious pubchempy. The source code was only restructured and safety data capabilities were added. Backward compatibility with the original. Also removed the python2 checks.

.. code:: python

    >>> from pubchempy import get_compounds, Compound, get_SDS
    >>> comp = Compound.from_cid(1423)
    >>> print(comp.isomeric_smiles)
    CCCCCCCNC1CCCC1CCCCCCC(=O)O
    >>> comps = get_compounds('Aspirin', 'name')
    >>> print(comps[0].safety_data)
    {'pictogram': [{'icon': 'GHS07.svg', 'string': 'Irritant'}, {'icon': 'GHS08.svg', 'string': 'Health Hazard'}], 'hazard': ['H302', 'H315', 'H316', 'H319', 'H334', 'H335', 'H360', 'H370', 'H371', 'H372', 'H373'], 'precautionary': ['P201', 'P202', 'P260', 'P261', 'P264', 'P270', 'P271', 'P280', 'P281', 'P285', 'P301+P312', 'P302+P352', 'P304+P340', 'P304+P341', 'P305+P351+P338', 'P307+P311', 'P308+P313', 'P309+P311', 'P312', 'P314', 'P321', 'P330', 'P332+P313', 'P337+P313', 'P342+P311', 'P362', 'P403+P233', 'P405', 'P501']}
    >>> print(pcp.request_SDS(cid=1254))
    {'pictogram': [{'icon': 'GHS07.svg', 'string': 'Irritant'}, {'icon': 'GHS05.svg', 'string': 'Corrosive'}], 'hazard': ['H315', 'H318', 'H319', 'H335', 'H402', 'H412'], 'precautionary': ['P261', 'P264', 'P271', 'P273', 'P280', 'P302+P352', 'P304+P340', 'P305+P351+P338', 'P310', 'P312', 'P321', 'P332+P313', 'P337+P313', 'P362', 'P403+P233', 'P405', 'P501']}





Installation
------------

Install PubChemPy using:

::

    git clone https://github.com/vinc3PO/PubChemPy.git
    cd to/path
    pip install .

Documentation
-------------

Full documentation of the original available at http://pubchempy.readthedocs.io.

Contribute
----------

-  Feature ideas and bug reports are welcome on the `Issue Tracker`_.
-  Fork the `source code`_ on GitHub, make changes and file a pull request.

License
-------

PubChemPy is licensed under the `MIT license`_.

.. _`installation options`: http://pubchempy.readthedocs.io/en/latest/guide/install.html
.. _`source code`: https://github.com/mcs07/PubChemPy
.. _`Issue Tracker`: https://github.com/mcs07/PubChemPy/issues
.. _`MIT license`: https://github.com/mcs07/PubChemPy/blob/master/LICENSE
