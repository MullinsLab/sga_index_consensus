## branch sga-index

This branch has been further adapted by Dylan Westfall to enable the pipeline to demultiplex SGA samples 
labeled with Index primers rather than unique Sample IDs present in the cDNA primer. 
This primarily revolved around replacing the demultiplexing code in demux_functions.jl with
code which Alec had written previously for an older version of PORPID-postproc that could accept
Indexed samples and removing code which handled UMIs and postprocessing reports. 


# PORPIDpipeline

by Alec Pankow and Ben Murrell, now maintained by Hugh Murrell

now upgraded to Julia version 1.7.1

## branch master 

## Quick start
## For installation on MacOS skip down to the Conda setup

### Dependencies (on an ubuntu machine)

- first update all apps
   - `apt update`
   - `apt upgrade`
- Snakemake
   - `apt-get install -y snakemake`
- mafft
   - `apt-get install -y mafft`
- fasttree
   - `apt-get install -y fasttree`
- python3 packages
  - `apt-get install python3-pandas`
  - `apt-get install python3-seaborn`


### Julia version 1.7

Download and unpack the latest Julia (we recommend version 1.7.1) from: 

[https://julialang.org/downloads/](https://julialang.org/downloads/)

Make sure you can enter the julia REPL from the command line, on an ubuntu machine you would do:

```bash
# move the julia system to a lib directory
mv julia-1.7.1 /usr/lib/julia-1.7.1
# make julia v1.7.1 executable from the command line
ln -s /usr/lib/julia-1.7.1/bin/julia /usr/local/bin/julia
# check that you can enter the julia REPL
julia --version
```

### cloning the PORPIDpipeline repository

Now that the dependencies are setup we clone the PORPIDpipeline repository

```bash
cd ~
git clone https://github.com/MurrellGroup/PORPIDpipeline.git
```

### setting up the Julia package environment

then navigate to the `PORPIDpipeline` project folder and start the Julia REPL. 
Enter the package manager using `]` and then enter

```julia
activate .
instantiate
precompile
```

This will activate, install, and precompile the `julia` environment specified by the 
`Project.toml` and `Manifest.toml` files. The `precompile` command
above is not strictly needed but is useful if there are issues with installing
the `julia` packages listed in `Project.toml`

If the following error occurs when using conda during instantiate,
"GitError(Code:ERROR, Class:SSL, Your Julia is built with a SSL/TLS
engine that libgit2 doesn't know how to configure to use a file or
directory of certificate authority roots, but your environment specifies
one via the JULIA_SSL_CA_ROOTS_PATH variable. If you believe your
system's root certificates are safe to use, you can `export
JULIA_SSL_CA_ROOTS_PATH=""` in your environment to use those instead.)"

leave the package manager and enter 

`ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""`

Return to the package manager `]` and try again. It is also possible for `instantiate` to fail
due to latency. If this happens rerun `instantiate` to try again.

Next, add the following text to your Julia startup file (typically at `~/.julia/config/startup.jl`; 
you may need to create the directory if not present, `mkdir -p ~/.julia/config`). For MacOS installations
with julia inside the PorpidPostproc environment, the config directory should be placed here instead.
`~/mambaforge/envs/PorpidPostproc/share/julia`

```julia
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
```

This will activate the local environment at Julia startup.


### Configuration

To configure the workflow, first edit the demo `config.yaml` file to reflect
your library construction. 
It should follow the same format as shown below in **config.yaml**

```yaml
Dataset1:
  Sample1:
    cDNA_primer: No cDNA primer listed for this sample
    fwd_index: Index_F21
    index_type: Index_primer
    panel: No panel listed for this sample
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample2:
    cDNA_primer: No cDNA primer listed for this sample
    fwd_index: Index_F22
    index_type: Index_primer
    panel: No panel listed for this sample
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample3:
    cDNA_primer: No cDNA primer listed for this sample
    fwd_index: Index_F23
    index_type: Index_primer
    panel: No panel listed for this sample
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
```

### Index Primers and Config Generation 
Each sequenced PCR product to be
used with this pipeline should have been indexed with one of the primer
combinations listed in Index_primer_sheet.csv file in the bash-scripts folder. The primer sequences provided will be used for demultiplexing and will be trimmed
from the final sequences. **sec_str_primer** and **rev_primer** are the PCR primers 
used in the 2nd rd of PCR. **fwd_index** and **rev_index** are the index primers used to 
tag each sample. Illumina Nextera indexes can be used as well as the Index primer set 
developed in Jim Mullins lab, which can be found in the bash-scripts folder.

To produce the config file required to run the pipeline for any given
sample set, a csv must first be created with the Index and PCR primer
information for each sample. It is critical to name this file correctly
as the name of this csv file is copied in as the dataset name in the
config file. This dataset name will need to match the CCS read file
(fastq.gz) and will be used to name the output data from the pipeline.

Example: 2024-01_FH13_bc1008.csv goes with dataset
2024-01_FH13_bc1008.fastq.gz

An example primer sheet can be found in the bash-scripts folder. (Dataset1.csv)

To create the config.yaml file for the sample set, the python script
primersheet2config_SGA-Index.py is run as described below.

primersheet2config_SGA-Index.py Located in the directory 'bash-scripts'.
Run this python script to generate a config.yaml file from a given csv.


Run: python primersheet2config_SGA-Index.py [dataset].csv >
[dataset]-config.yaml

(Example:  python primersheet2config_SGA-Index.py Dataset1.csv >
Dataset1-config.yaml)

After creating the config file it must be moved to the main sga-index-consensus
directory or the contents copied into the existing config.yaml file.

This config.yaml file is the main config file used per
pipeline run, dictating the dataset name, what samples are included, and
also any primer sequences and panel files required for each sample. This
is created from the primersheet2config_SGA-Index.py file, as described
above, and moved to the sga-index-consensus folder or copied to the
existing config.yaml file. Double check this config file looks correct
before a run to ensure primers are correct. A sample which gives an
issue or crashes the pipeline can be commented out to prevent the
pipeline from attempting to run it.

The snakefile contains the rules for the pipeline and also
parameters that may need to be changed depending on the run. For
example, with recent NFLG samples, the max_length parameter was changed
to 10,000 (from the default 4000) and min_length changed to 10 (from the
default 2100). This file mostly remains the same, other than any
specific snakefile parameter changes. Be sure that the snakefile is
pointing to the correct config file (line 3, `configfile:
"config.yaml"`) or that the correct contents have been copied to the
config.yaml file.

Raw CCS .fastq files should be placed in the `raw-reads/` subdirectory and named 
according to the the dataset name used in the `config.yaml` file, ie, `Dataset1.fastq`
for the example dataset. For example, a config file with dataset
name "2024-01_FH13_bc1008" will only work with a raw read file named
"2024-01_FH13_bc1008.fastq.gz".

The contamination check will compare all sequences to one another, but if it is desired for
another sequence (for example a lab strain) to be included any sequence can be added to the
`contam_panel.fasta` file in the panel folder. The sequence must be trimmed to the same length
as the target sequence it is to be compared to. Note that the contamination check will flag
sequences that are within the distance threshold `dist_thresh` specified in the snakefile
(default 1.5%).

### Pipeline Results
The output files from the pipeline are all written to dataset/[dataset
name]. There are 4 summary reports and 3 folders containing output from
the pipeline. The reports are described below:

reject_report.csv - the number of reads rejected by the filter for CCS
read quality 

quality_report.csv - a table summarizing total reads,
quality reads, and total number of reads assigned to samples

demux_report.csv- the number of reads matching the index and PCR primers
for each sample in the config. 

contam_report.csv - all consensus sequences
are compared to one another and any that are <1.5% different from
another sequence are added to the table and the non-self sample and
distance are indicated in the table. Be aware that sequences from the
sample individual will nearly always be flagged due to similarity. 

The 3 folders containing sequence outputs from the pipeline are summarized
below:

Demux - fastq files containing reads with all index and PCR primers
trimmed off. Samples with no reads found will not have a file here
(consult demux_report.csv). Any sample with reads where the PCR primer
sequences do not match the config are written to files with the suffix
"_primer_trim_issue.fastq". These have not had any sequence trimmed off
so Index and PCR primer can be checked manually.
REJECTS_DEMUX_QUAL.fastq contains reads which failed the quality filter
before demultiplexing began. 

Consensus - consensus sequences (.fasta) created from the reads in each
fastq file. Sequence names have additional information added to the
names: 
 1) the number of reads used to create the consensus (num_CCS=) 
 2) the minimum agreement at any position inthe alignment of reads used to create the consensus sequence
   (min_agreement=) 
   
Filtered - consensus sequence files (.fasta) are filtered into two
folders by whether their minimum agreement value is above or below 0.7.
Files in the passed_0.7_agreement folder are all consensus sequences
with min_agreement > 0.7. Any sequence with min_agreement > 0.7 is
considered to be high quality and analyses ready. Sequences in the
folder "below_0.7_agreement" all have consensus sequences with
min_agreement > 0.7. Any sequence with min_agreement < 0.7 and need to
be reviewed manually to determine the cause of the low agreement. To aid
this effort, alignments for each sample below the minimum agreement are
written to the alignments folder. 

Common causes of minimum agreement failure are listed below. Depending on
analysis requirements some or all of these may be acceptable. By
consulting the alignment a cause can be determined and the sequence
confirmed unusable or perhaps edited as appropriate and included in
downstream analysis. 

1) long homopolymer runs where PCR or sequencing polymerases slip and
create variable lengths among reads,

2) 1st cycle PCR errors, where an error is made by the polymerase when
copying cDNA, leading to 50/50 mix of nucleotides at that position in
amplicons derived from that double stranded DNA

3) Multiple templates in the PCR well, leading to differing frequencies
at any position of disagreement


### Preview and Execution

Preview jobs with Snakemake and run with {n} cores.

```bash
#preview jobs
snakemake -np

#run
snakemake -j{n}
```

For more info on Snakemake, see https://snakemake.readthedocs.io/en/stable/

## Conda setup

Some (without root access) may prefer to setup PORPIDpipeline in a **conda** environment.

To accomplish this, first install `anaconda` locally. (the install script allows you to choose
the location for anaconda, by default `/home/user` but choose something else if
you want something accessable to a group of users)

```bash
curl –O https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh > Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```

Or for MacOS visit the anaconda documentation and install the most recent version
https://docs.anaconda.com/anaconda/install/mac-os/

then log out and log in again and check that you are in the `base` environment.

`conda` is very slow, so we suggest installing `mamba` in the conda `base` environment:

Visit the mamba documentation and download the most recent version
https://mamba.readthedocs.io/en/latest/installation.html
https://github.com/conda-forge/miniforge#mambaforge

For MacOS
```bash
cd Downloads
bash Mambaforge-MacOSX-x86_64.sh
```

When asked if you "wish the installer to initialize Mambaforge
by running conda init? [yes|no]" 

select yes

then log out and log in again and check that you are in the `base` environment.

clone the PORPIDpipeline repository

```bash
cd ~  # or some other directory used for your anaconda installation
git clone https://github.com/MurrellGroup/PORPIDpipeline.git
```

and then all the PORPIDpipeline dependencies including `julia` version `1.7.1`
( as listed in the PORPIDpipeline conda environment spec in `environment.yaml`),
can be installed in a `conda` environment via `mamba` using the commands:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
mamba env create --file environment.yaml
```

Note that if you did use *some other directory* than your home directory for
installing the PORPIDpipeline repository then you have to inform Julia where
your packages are stored by placing the following command in your `.bashrc`
file. If installing on MacOS skip this step:

```bash
# set path to .julia files
export JULIA_DEPOT_PATH="/some/other/directory/.julia"
```

to complete the setup, activate the new PORPIDpipeline conda environment, 
#currently environment.yaml title is PorpidPostproc instead of PORPIDpipeline

```bash
conda activate PorpidPostproc
```

and continue with the `julia` package environment setup as outlined above in the *quick start* section.

## Cluster setup

Seting up a `snakemake` pipeline on a cluster is a *dark* art. Here we describe an attempt
at installing PorpidPostproc on a two node cluster, (one node a *controller* node with 16 cores
and the other node a *compute* node with 64 cores).

**Firstly**, since the cluster administrator is hardly likely to give you root access we
suggest you follow the `conda` installation for PorpidPostproc. If you expect more
than one user of your PorpidPostproc pipeline then install in a directory that all
your users can see and that is visible from both the *contoller* and *compute* nodes. 
ie use `some other directory` rather than the standard home directory and make 
sure to inform `julia` about this choice of directory as
outlined in the `conda` section above.

**Secondly**, cluster administrators usually insist that large data sets are stored
in an appropriate volume and **not** in the usual user's space. On our cluster the
administrator required the PorpidPostproc code to be installed in a `\tools\porpid\`
directory and the large data sets (input, output and temporary) to be stored in
a `\data\porpid\` directory so we installed PorpidPostproc into `\tools\porpid\porpidpostproc`
and then replaced some of the directories in the `porpidpostproc`
directory with symbolic links to an appropriate directory in the `\data\porpid\` directory
as shown below

```
config.yaml -> /raw/porpid/config/demo.yaml
panels -> /raw/porpid/panels/
porpid -> /raw/porpid/porpid/
postproc -> /raw/porpid/postproc/
raw-reads -> /raw/porpid/raw-reads/
```

Naturally, one must copy contents of the installation to the `/raw/porpid/` directory
before deleting the installation directory and replacing it with a symbolic link to the
appropriate place on the `raw` volume.

**Job submission**, after setting up like this we are ready to run the `demo` study through PorpidPostproc
by submitting the `snakemake` command to the cluster managemant system.
On our cluster that management system is `slurm` and the following shell script
stored in `porpid_job.sh` facilitated that submission:

```bash
#!/bin/bash
#SBATCH --job-name==porpid
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --partition=main

if [ "$#" -lt 1 ]; then
    echo "please supply a config file name as first parameter"
    exit
fi
echo "config file is $1"

echo "${SLURM_JOB_NAME} job submited using ${SLURM_NTASKS} cores"

# create a symbolic link for the snakemake config file to point to the config for the current study
rm -f /tools/PorpidPostproc/porpidpostproc/config.yaml
ln -s /RAW/PORPID/CONFIG/$1.yaml /tools/PorpidPostproc/porpidpostproc/config.yaml

# tell slurm where anaconda is and conda activate the PorpidPostproc environment
source /tools/PorpidPostproc/anaconda3/etc/profile.d/conda.sh
conda activate PorpidPostproc

# navigate to the porpidpostproc directory and run snakemake
# add -F to to the snakemake command to force re-run of all rules
cd /tools/PorpidPostproc/porpidpostproc
snakemake --rerun-incomplete -j${SLURM_NTASKS}  
```

To submit the `demo` to run as a `slurm` batch job one just uses

```bash
sbatch porpid_job.sh demo
```
The script above sets some environment variables for `slurm` and then resets
the symbolic link to the appropriate config file for the `demo` study.
It then activates the conda environment switches to the installation
directory and runs the snakemake pipeline.

With this structure it is easy to run a new study through PorpidPostproc.
One copies the new config file into the `/raw/porpid/config/` directory,
transfers the `fastq` data to the `/raw/porpid/raw-reads/` directory
and then issues the `sbatch` command using the appropriate study name
instead of `demo`

Note that with this method you must predetermine the number of cores
you intend to use on your cluster's node. In the `demo` study this is set
to 7 ( 6 cores for the samples to run in parallel plus 1 core for snakemake )

Each study will be different. To see how many samples can be run in parallel
you can do a `snakemake` dry run using the `porpid_dry_run.sh` script below:

```bash
#!/bin/bash
if [ "$#" -lt 1 ]; then
        echo "please supply a config file name as first parameter"
        exit
fi
echo "config file is $1"
# create a symbolic link for the snakemake config file to
# point to the config for the current study
rm -f /tools/PorpidPostproc/porpidpostproc/config.yaml
ln -s /RAW/PORPID/CONFIG/$1.yaml /tools/PorpidPostproc/porpidpostproc/config.yaml
# activate the conda environment
source /tools/PorpidPostproc/anaconda3/etc/profile.d/conda.sh
conda activate PorpidPostproc
# perform a snakemake dry run
# remove the -f for a partial dry run of what's left to do
cd /tools/PorpidPostproc/porpidpostproc
snakemake -F --rerun-incomplete -np
```

Note that this dry run is not compute intensive and can ve executed on the
*controller* machine without using the `sbatch` command as follows:

```bash
./porpid_dry_run.sh demo
```

### Caveat

The above suggestion for running a `snakemake` pipeline under `slurm`
is rudamentary. Maximum cores must be requested at the start of execution
and they are probably held throughout the run.

However, it is alledged that `snakemake` can play nicely with `slurm` and
it should be possible to have `snakemake` invoke `slurm` for each rule in
the pipeline. In this case `snakemake` would request the optimal number
of cores needed for each step in the pipeline.

We have not attempted this yet, and it would probably require writing a
`slurm` efficient version of the `snakefile`. 

Watch this space for further developments.

## Documentation

### Workflow

The graph below summarizes the overall organization of the workflow. 
Each node in the graph is a *rule* in the The [Snakefile](Snakefile).

![rulegraph](rulegraph.png)





