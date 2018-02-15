# CRISPR_immunity
Data and code for paper arXiv:1801.10086


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
