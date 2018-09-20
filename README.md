# PhyloTempo
Temporal clustering of a phylogeny

**This project is still in beta and as such may contain major bugs and we do not guarantee proper results/functionality**

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
temporalClustering("file.tree", "timetable.txt")
```

Important functions (with defaults):

```R
temporalClustering(tree_file, timetable_file, parsimony=FALSE , bootstrap=200 , output="TC" , randomMulti2Di=TRUE )
```

## Output

PhyloTempo generates two figures in PDF format and one table:

<object data="./img/Ancestral_character_tree.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="./img/Ancestral_character_tree.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="./img/Ancestral_character_tree.pdf">Download PDF</a>.</p>
    </embed>
</object>

<object data="./img/TC_plots.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="./img/TC_plots.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="./img/TC_plots.pdf">Download PDF</a>.</p>
    </embed>
</object>

 | Set | TimeRange | TimeIntervals | Tips | PearsonRho | TC | StaircaseNess | Cherries | PybusGamma | Colless | Sackin | CollessYule | CollessUnif | SackinYule | SackinUnif | 
 | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
 | TC | 1 - 9 time_units | 9 | 137 | 0.864449781427698 | 0.247069329881309 | 0.661764705882353 | 46 | -35.5760189899849 | 557 | 961 | 0.261644025388755 | 0.347355631641048 | -1.98708503154303 | 0.599297597858254 | 

## References
Norström MM, Prosperi MC, Gray RR, Karlsson AC, Salemi M (2012). Phylotempo: a set of r scripts for assessing and visualizing temporal clustering in genealogies inferred from serially sampled viral sequences. Evol Bioinform Online. 8:261.
