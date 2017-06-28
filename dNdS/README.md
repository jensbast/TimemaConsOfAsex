## Introduction

dN/dS ratio estimates of Timema using [Codeml (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html).

Alignments are in [phylip](http://evolution.genetics.washington.edu/phylip.html) format.

Trees have been generated with [RAxML](http://sco.h-its.org/exelixis/software.html) but have to be reformatted depending on the model.

```timema2codeml.py``` is useful to reformat the input data to make it compatible with [Codeml (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html).
 
## Table of contents

1. [Required files](#1_input)
2. [dN/dS ratio estimates](#2_paml)

## <a name="1_input"></a>1) Required input

1. Need 2 directories and 1 control file:

	* alignments/ directory: contain alignments in phylip format.
	* branch_lengths/ directory: contain the RAxML's trees in newick format.
	* control template (codeml_timema.ctl)

## <a name="2_paml"></a>2) dN/dS ratio estimates

1) Create a directory for your specie (timema).

```
mkdir timema
```

Copy the 3 required input into this directory.

```
ll timema

drwxr-xr-x 2 ptranvan dee_schwan 159744  2 mar 11:00 alignments
drwxr-xr-x 2 ptranvan dee_schwan 241664  2 mar 11:01 branch_lengths
-rw-r--r-- 1 ptranvan dee_schwan 3070  2 mar 11:00 codeml_timema.ctl
```

```timema2codeml.py``` can be stored anywhere but has to be run from the ```timema``` directory.

2) Run the script for 3-rate, 2-rate and free models.

```
cd timema

python timema2codeml.py -s timema -i1 alignments -i2 branch_lengths -i3 codeml_timema.ctl -e <your_email> -o paml.sh
```

3) A script (*.sh) will be created. It has to be run on Vital-IT server (usr@prd.vital-it.ch):

```
chmod +x *.sh
./*.sh
```

Once done, a directory output will be created.

4) You can create a table for each model to summarize the results.

```
python timema2codeml.py -s table_3rate -i1 output/ -o 3rate_table.txt

python timema2codeml.py -s table_2rate -i1 output/ -o 2rate_table.txt

python timema2codeml.py -s table_free -i1 output/ -o free_table.txt
```

Annexe) Options:

```
-i1: alignments directory.
-i2: branch_lengths directory.
-i3: control file.
-e: email for job confirmation.
-o: name of the script file to be executed.
```

