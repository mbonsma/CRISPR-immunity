# CRISPR immunity

Data and code for paper [arXiv:1801.10086](https://arxiv.org/abs/1801.10086): **How adaptive immunity constrains composition and fate of large bacterial populations**

### Data

We analyzed data collected by Paez-Espino *et al.* (2013) and reported in [this paper](https://www.ncbi.nlm.nih.gov/pubmed/23385575). Original data from their study is publicly available in the NCBI Sequence Read Archive under the accession [SRA062737](https://www.ncbi.nlm.nih.gov/sra/?term=SRA062737). It includes four data files (SRR630110, SRR630111, SRR630412, and SRR630413) which we used for our analysis. We extracted the data corresponding to the MOI2 deep sequencing experiment and separated it into time points by checking each read for matches to the primers identified in the supplementary information of Paez-Espino *et al.* 

The data files deposited here are a result of our further processing and sorting.

 * [CRISPR1_MOI2_deep_spacer_types.txt](https://github.com/mbonsma/CRISPR_immunity/blob/master/data/CRISPR1_MOI2_deep_spacer_types.txt) contains every unique spacer sequence we detected and its assigned type. Spacers were grouped such that sequences within an edit distance of 2 were assigned to the same type.
	* column 1: spacer sequence
	* column 2: spacer type label

 * [CRISPR1_MOI2_deep_timepoint_1.txt](https://github.com/mbonsma/CRISPR_immunity/blob/master/data/CRISPR1_MOI2_deep_timepoint_1.txt) contains spacers arranged by their position in each read. Each row represents a read, and the columns contain spacers sequences if one was present or zeros if no spacers are present in that position. Positions are ordered such that position 0 is next to the wildtype spacers and position 5 is the furthest from wildtype (newest) spacer. The remainder of the files in the `data` folder are the same structure as this file but for the rest of the time points in the study (2 through 14). 

**NOTE**: time point 8 was nearly completely missing in the original data, which was confirmed through communication with the authors. We excluded time point 8 from our analysis.  

We also analyzed data from [Parikka *et al.* (2017)](https://www.ncbi.nlm.nih.gov/pubmed/27113012) to compare virus-to-prokaryote ratios (VPR) from a variety of environments with our model. Their data is deposited at [10.15454/1.4539792655245962E12](https://www.researchgate.net/publication/312027517_Data_and_metadata_dealing_with_prokaryote_and_viral_abundances_from_a_variety_of_ecosystems)

### Compiling and running simulation code

Simulation code written by Dominique Soutiere.

 1. Install the [`Eigen`](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) package.

 2. To compile the code:

``` 
g++ -std=c++0x -I ~ -O3 tauLeapingWithLoss3C0.cpp -o tauLeaping
```

After the `-I` flag, put the path to the `Eigen` folder downloaded in step 1. If it's in the home directory, use `-I ~`, otherwise use `-I path/to/Eigen`. 

 3. Make a `./data/` folder in the same folder as the compiled C++ code to save the output data to.

 4. To run the code: 

```
./tauLeaping m xi eta R pv0 exp c0Exp repeat initial_x initial_y initial_z initial_nu
```

Input parameters:

 * `m` total number of possible spacer types
 * `xi` = $1-e$, where $e$ is the spacer effectiveness ranging from $0$ to $1$. 
 * `eta` spacer acquisition probability ($\leq 1$)
 * `R` = $r/(gC_0)$, spacer loss rate 
 * `pv0` = $p_V$, probability of phage success against bacteria without spacers
 * `exp` maximum iterations for simulation = $10^{\text{exp}}$
 * `c0Exp` $C_0$, nutrient concentration = $10^{\text{c0Exp}}$
 * `repeat` replicate indicator, used only in creating output filename
 * `intial_x` initial bacterial population size ($0 \leq x \leq 1$), $x = n_B/C_0$
 * `initial_y` initial phage population size ($y \geq 0$), $y = n_V/C_0$ 
 * `initial_z` initial nutrient concentration ($0 \leq z \leq 1$), $z = C/C_0$
 * `initial_nu` initial fraction of bacteria with spacers ($0 \leq \nu \leq 1$), $\nu = n_B^s/n_B$


Example usage:

```
./tauLeaping 500 0.5 0.0001 0.04 0.02 5 7 0 0.3 3 0.7 0
```
