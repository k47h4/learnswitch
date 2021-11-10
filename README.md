requires python 3.6, pypet, matplotlib, seaborn and MATLAB_R2020A

install miniconda: https://docs.conda.io/en/latest/miniconda.html
create a conda environment with python version 3.6:
>conda create --name <environment name> python=3.6
activate this environment:
>conda activate <environment name>
install the required modules:
>conda install pip
>pip install pypet
>conda install matplotlib
>conda install seaborn


To reproduce the data in Fig. 6B and C for the multiplicative model:
run the file as is:
>python run_model.py
(The parameter 'modulation' in line 459 should be 'multiplicative'
the variable 'identifier' in line 601 should be 'multiplicative')

To reproduce the data in Fig. 6C for the additive model:
Modify the file run_model.py:
Change the parameter 'modulation' in line 459 to 'additive'
Change the variable 'identifier' in line 601 to 'additive'

To reproduce the data in Fig. 6D and E:
Modify the file run_model.py:
Change the parameter 'modulation' in line 459 to 'multiplicative'
Change the variable 'identifier' in line 601 to 'pyrsst'

Once the hdf5 files are generated, run the following scripts to generate the mat files:

Figure 6B:
>python plot_responses.py

Figure 6C:
>python plot_matrices.py

Figure 6E and Figure S3:
>python plot_Fig6DE.py


To get the original figures from the paper, run the Matlab file in MATLAB:
plot_Figure6.m


The scripts take the data from the hdf5 files and mat files in the local folder and plot the figures into the local folder.

