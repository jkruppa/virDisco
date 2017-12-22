# Changes to the Pauda source code

Several changes to the Pauda implementation were necessary
to fulfill our requirements. The following Table shows all
changes with the connected program line. First, we changed the command
`-Xmx`, which specified the maximum memory allocation pool to 50 GB
in all shell scripts. This was needed to allow Pauda to handle the large
amino acid reference genome. Second, we changed some Bowtie2 parameters. We
set `-N`, the number of mismatches to be allowed in a seed alignment
from 1 back to the default value of 0. This makes the alignment less
sensitive but much faster. We set `-K` from 30 to 100, allowing
Bowtie2 up to one hundred multi maps per read. Finally, we set
`-score-min` from 35 to 10, which is the minimum score for a
alignment to be judged as valid.

| Source file  | Line          | From  | To   |
|------------- |:-------------:| -----:|-----:|
| z1\_protein2pna.sh | 56 | --Xmx20 | --Xmx50|
| z4\_bowtie-on-pna.sh | 26 | bowtie2-align | bowtie2|
| | 45 | N=1 | N=0|
| | 46 | K=30 | K=100|
| | 49 | MS=35 | MS=10| 
| z5\_sam2blastx.sh | 50 | --Xmx8 | --Xmx50|
| | 59 | --Xmx8 | --Xmx50|
| z6\_blastx2rma.sh | 27 | --Xmx8 | --Xmx50|

