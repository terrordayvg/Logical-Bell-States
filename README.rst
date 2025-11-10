=========================================================================================================================
``Feasibility of Logical Bell State Generation in Memory Assisted Quantum Networks``
=========================================================================================================================

.. image:: https://img.shields.io/badge/python-3.11-blue.svg
        :target: https://www.python.org/downloads/release/python-3110/

.. image:: https://img.shields.io/badge/arXiv-2412.01434-b31b1b.svg
        :target: https://arxiv.org/abs/2412.01434


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
We conduct a feasibility analysis of employing logical Bell states in quantum memories for quantum networks. Specifically, we determine the break-even and pseudo-thresholds for Surface and Bacon-Shor codes, focusing on codes of distance \(d=3\) and \(d=5\). To generate the logical Bell states, we utilize lattice surgery and introduce two memory assisted protocols: one employing local generation and the other utilizing non-local generation.


How to Cite
===========

If you use this work, please cite the following papers:
::

    @article{rsrk-c7yg,
       title={Feasibility of logical Bell state generation in memory assisted quantum networks},
       author={Galetsky, Vladlen and Vyas, Nilesh and Comin, Alberto and Nötzel, Janis},
       journal={Phys. Rev. Res.},
       volume={7},
       issue={3},
       pages={033090},
       numpages={14},
       year={2025},
       month={Jul},
       publisher={American Physical Society},
       doi={10.1103/rsrk-c7yg},
       url={https://link.aps.org/doi/10.1103/rsrk-c7yg}
    }
