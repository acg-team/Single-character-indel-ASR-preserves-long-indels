# The Python scripts used for the paper entitled:
### Single-character insertion-deletion model preserves long indels in ancestral sequence reconstruction

## Requirments
For this tutorial, you should already have [Python 3.9 or higher](https://realpython.com/installing-python/), [jupyter notebook](https://jupyter.org/install) along with the following libraries:

```numpy, ete3, biopython, seaborn, matplotlib, sklearn and pandas```

### Installation

To install the package you can simply download the repository and run the following command in the root directory.

Install the dependencies using this command:

```console
pip3 install -r requirements.txt
```

### Files

<ul>
	<li> mammals_01_stat.ipynb contains functions for indel pattern plots for mammalian data.</li>
	<li> mammals_02_dynamic_of_gaps.ipynb includes functions for calculating dynamic of gap pattern for each mammalian data.</li>
	<li> mammals_03_indel_length.ipynb includes functions for ploting indel length for mammalian data.</li>
	<li> simulation_01_acc.ipynb contains functions for computing accuracy of ARPIP inference on simulated data.</li>
	<li> simulation_02_dynamic_of_gaps.ipynb contains functions for calculating dynamic of gap pattern for each simulated data. </li>
	<li> simulation_03_stat.ipynb contains functions for indel pattern plots for simulated data.</li>
	<li> simulation_04_discussion.ipynb contains scripts for the appendix figures.</li>
	<li> requirements.txt contains library versions of dependencies.</li>

</ul>
To get the figures in the manuscript all the necessary files and scripts are provided here.
Moreover, suplemental data is stored in another repository with [this link](https://doi.org/10.5281/zenodo.10798097). 

## Citation

Please cite:

Gholamhossein Jowkar, Julija Pecerska, Manuel Gil, and Maria Anisimova  
**Single-character insertion-deletion model preserves long indels in ancestral sequence reconstruction.**  
BioRxiv, 2024;  
doi:[10.1101/2024.03.09.584071](https://doi.org/10.1101/2024.03.09.584071)

---
## Author

Gholam-Hossein Jowkar [E-mail](jowk@zhaw.ch)