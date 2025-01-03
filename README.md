SyntenyLink 🧬
=============

## Table of Contents 📚

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Main programs](#main programs)
6. [Contact](#contact)
7. [License](#license)

## Overview 📖
===========

The SyntenyLink package has six major components: the SyntenyLink
algorithm allows users to handle reconstruct subgenomes of polyploid
species more conveniently and to separate the set of genes belong to
each subgenome in the organism with the aid of reference proteomes of
polyploid species and related ancestor. 🌱

For more details access our published work here: https://ieeexplore.ieee.org/document/10385622

![Pipeline](figs/pipeline_Algo.png)

All programs are executed using command line options on Linux systems or
Mac OS. Usage or help information are well built into the programs. 💻

All code is copiable, distributable, modifiable, and usable without any
restrictions. 

## Requirements 🛠️
=============

To use SyntenyLink, ensure you have the following requirements and python packages:
### Python packages:

    - Python
    - biopython
    - ipython
    - matplotlib
    - numpy
    - pandas
    - seaborn
    - pickle
    - csv
    - os
    - math
    - sys
    - re
    - warnings

### Other tools:

    - makeblastdb
    - blastall
    - dagchainer

## Installation ⚙️
=============

1. Clone this repository to your local machine:

```bash
git clone https://git.cs.usask.ca/qnm481/syntenylink.git
cd syntenylink
```

2. Install the required dependencies:

```bash
pip install -r requirements.txt
```

## How to use SyntenyLink 🚀
=============

3. Reproduce all the experiments:

i. Run SyntenyLink.sh
```
./SyntenyLink.sh ref_pep.fasta ref_cds.fasta query_pep.fasta query_cds.fasta query.gff3 ref.gff3 ref_genelist.txt query.bed -n <number of subgenomes> -s <ploidy status> -chr1 <query chromosome number for subgenome1> -p <gene prefix> 
```

## More information 🚀
=============

Utilize this repository to replicate our experiments and explore the functionalities of SyntenyLink. The codebase is organized to help you easily navigate through different components and reproduce our results.


The following is the list of executable programs
------------------------------------------------

## Usage
=============

Parameter and command examples are shown below.

1. [SyntenyLink.sh](#syntenyLink.sh) 
    - -i Input collinear file
    - -s Ploidy status. If diploid 2, tetraploid 4, hexaploid 6, octaploid 8
    - -p Gene prefix. 
    - -gt Groundtruth subgenome separation. Default is None
    - -bed Query bed file
    - -chr<i> Chromosome number for each subgenome. For example in Brassica napus chr1 is 10 and chr2 is 9
    - -n number of subgenomes.
    - -sub<i> Prefix for subgenome chromosome name. For example in Brassica napus sub1 is A and sub2 is C
    - -dag Collinear file created with blastn output

## Tested species
=============

1. _Brassica rapa_
2. _Brassica oleracea_
3. _Brassica nigra_
4. _Brassica napus_
5. _Brassica carinata_
6. _Brassica juncea_
7. _Sinapis alba_
8. _Cercis canadensis_

## Contact 📬
===============

For any questions or inquiries, please feel free to open an issue on our repository or contact us at [qnm481@usask.ca](mailto:qnm481@usask.ca).

## License 📜
===============

This project is licensed under the [MIT License](LICENSE)

