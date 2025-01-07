=========================================================================================================================
``Feasibility of Logical Bell State Generation in Memory Assisted Quantum Networks`` 	arXiv:2412.01434
=========================================================================================================================


.. image:: https://codecov.io/gh/terrordayvg/PT_VQC-Tomography/graph/badge.svg?token=880RTY0T96
        :target: https://codecov.io/gh/terrordayvg/PT_VQC-Tomography

.. image:: https://img.shields.io/badge/python-3.11-blue.svg
        :target: https://www.python.org/downloads/release/python-3110/



Installation of required libraries

::

    install -r requirements.txt


Usage

               This repository is divided into two main folders of the static and robust noise models
        
        Folders and files:  
                * Each folder has the images to be recreated in the article.
                * Has a sim.py file to produce the dataset or images.
                * And a plot.py file to recreate the exact figures with the dataset used.

                
        Aditional: 
                * IBM Sherbrooke 127-qubit calibration data is necessary and is presented in a .csv file.


Contents of requirements.txt
-----

::     

        gen==0.1
        matplotlib==3.5.2
        numpy==2.2.1
        PyMatching==2.1.0
        pytest==7.4.2
        qiskit==1.2.4
        qiskit_aer==0.15.1
        qiskit_ibmq_provider==0.19.2
        qiskit_ignis==0.7.1
        qiskit_terra==0.25.2.1
        scipy==1.15.0
        sinter==1.12.1
        stim==1.12.1
        qiskit_terra==0.22.3



What we aimed with this work?
To create a feasibility analysis of using logical bell states in quantum memories for quantum networks. We determine the tresholds and pseudo-tresholds of the code and code families for Surface and Bacon-Shor codes (d=3 and d=5 only). We utilize lattice surgery to generate the logical bell states using two different protocol with a local and a non-local generation.

-------------------

This repository combines all the codes to produce the plots and results from the following article: arXiv:2412.01434, if used, cite it correspondently. 
