# PhyloTempo (Beta)
Temporal clustering of a phylogeny

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

## Summary

This script calculates a number of [tree statistics](https://github.com/ProsperiLab/PhyloTempo#table) including the *temporal clustering* (TC) statistic which assesses the temporal structure of a phylogenetic tree by evaluating the order of time changes from internal nodes to tips.

## Requirements
R 3.3 or higher

**Libraries**

* ape
* doBy
* phytools
* quantreg
* diversitree
* apTreeshape
* RColorBrewer

The required input of PhyloTempo is a phylogenetic tree file in “newick” format and a two-column text file in which each tip name present in the phylo-genetic tree is associated with its corresponding time of sampling (a numeric value such as days or years).

## Installation

While most packages are easily installed:

```R
install.packages(c("ape","doBy","phytools","quantreg","apTreeShape","RColorBrewer"))
```

The package diversitree may give you an error and ask if GSL is installed. 
If it does, this may be easily resolved.

On Ubuntu this may be solved by first installing GSL:

```bash
sudo apt-get install libgsl0ldbl
```

For Centos/RHEL:

```bash
yum install gsl-devel
```

For Mac, you may need to install MacPorts, Homebrew or Fink and then install GSL.

```bash
brew install gsl # Homebrew
port install gsl # MacPorts
fink install gsl # Fink
```

## Execution 

```R
source("temporalClustering.R")
temporalClustering("Shankarappa.tree", "Shankarappa.txt")
```

Main function (with defaults):

```R
temporalClustering(tree_file, timetable_file, parsimony=FALSE , bootstrap=200 , output="TC" , randomMulti2Di=TRUE )
```

The Shankarappa files were originally produced in Norström *et al.* 2012 (see [references](https://github.com/ProsperiLab/PhyloTempo#references)).

## Output

PhyloTempo generates two figures in PDF format and one table:

### PDFs

<object data="./img/Ancestral_character_tree.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="./img/Ancestral_character_tree.pdf">
        <p>View <a href="./img/Ancestral_character_tree.pdf">Ancestral Character Tree PDF</a>.</p>
    </embed>
</object>

<object data="./img/TC_plots.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="./img/TC_plots.pdf">
        <p>View <a href="./img/TC_plots.pdf">Temporal Clustering PDF</a>.</p>
    </embed>
</object>

### Table
The table contains 15 columns:

- **Set** - Output name provided in method
- **TimeRange** - range of time units used for temporal clustering
- **TimeIntervals** - number of time intervals used for temporal clustering
- **Tips** - tip count
- **PearsonRho** - the Pearson correlation coefficient
- **TC** - The TC statistic assesses the temporal structure of a phylogenetic tree by evaluating the order of time changes from internal nodes to tips
- **StaircaseNess** - the proportion of sub-trees that are imbal-anced (ie, sub-trees where the left child contains more leaves than the right child, or vice-versa) compared against the distribution of such proportions obtained from random trees
- **Cherries** - the number of pairs of leaves that are adjacent to a common ancestor node
- **PybusGamma** - result of γ-test of Pybus & Harvey (2000)
- **Colless** - the Colless imbalance number of the tree
- **Sackin** - the Sackin index of the tree
- **CollessYule** - Colless' shape statistic under the Yule hypothesis
- **CollessUnif** - Colless' shape statistic under the Uniform hypothesis
- **SackinYule** - Sackin's shape statistic under the Yule hypothesis
- **SackinUnif** - Sackin's shape statistic under the Uniform hypothesis

 | Set | TimeRange | TimeIntervals | Tips | PearsonRho | TC | StaircaseNess | Cherries | PybusGamma | Colless | Sackin | CollessYule | CollessUnif | SackinYule | SackinUnif | 
 | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
 | TC | 1 - 9 time_units | 9 | 137 | 0.864449781427698 | 0.247069329881309 | 0.661764705882353 | 46 | -35.5760189899849 | 557 | 961 | 0.261644025388755 | 0.347355631641048 | -1.98708503154303 | 0.599297597858254 | 

## References
Bortolussi N, Durand E, Blum M, & François O (2006). apTreeshape: statistical analysis of phylogenetic tree shape. Bioinformatics. 22:3 https://doi.org/10.1093/bioinformatics/bti798

FitzJohn, RG (2012). Diversitree: Comparative Phylogenetic Analyses of Diversification in R. Methods in Ecology and Evolution (in press). doi:10.1111/j.2041-210X.2012.00234.x

Norström MM, Prosperi MC, Gray RR, Karlsson AC, & Salemi M (2012). Phylotempo: a set of R scripts for assessing and visualizing temporal clustering in genealogies inferred from serially sampled viral sequences. Evol Bioinform Online. 8:261.

Paradis E, Claude J & Strimmer K (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.

Revell, LJ (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223. doi:10.1111/j.2041-210X.2011.00169.x
