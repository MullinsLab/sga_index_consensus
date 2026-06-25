## branch sga-pipeline

This branch has been adapted by Dylan Westfall to enable the pipeline to
demultiplex SGA samples labeled with Index primers rather than unique
Sample IDs present in the cDNA primer. This primarily revolved around
replacing the demultiplexing code in demux_functions.jl with code which
Alec Pankow had written previously for an older version of
PORPIDpipeline that could accept Indexed samples and removing code which
handled UMIs and postprocessing reports.

Based upon code from PORPIDpipeline by Alec Pankow and Ben Murrell, now
maintained by Hugh Murrell.
(https://github.com/MurrellGroup/PORPIDpipeline.git)

now upgraded to Julia version 1.10.5

## Quick start


### Dependencies 

For an ubuntu machine:
- first update all apps
   - `apt update`
   - `apt upgrade`
- Snakemake and mafft
   - `apt-get install -y snakemake`
   - `apt-get install -y mafft`
- python3 packages
  - `apt-get install python3-pandas`
  - `apt-get install python3-seaborn`
  
### Julia version 1.10.5

We recommend you use the `juliaup` version manager to install julia.
from a terminal you can do this as follows:

```bash
curl -fsSL https://install.julialang.org | sh
```

This should install the Julia version manager, `juliaup` as well as
the latest version of Julia. To find out how to use the version manager 
to makesure you have version 1.10.5 as your default, go here:

[https://github.com/JuliaLang/juliaup]

Once Julia is installed, make sure you can enter the julia REPL from 
the command line and check the version number by logging out and in again 

```bash
exit
ssh root@.....
```

and then from your new terminal session:

```bash
juliaup status
julia --version
```

If the version number is not 1.10.5 then you need to use `juliaup` to install
that version and make it the default. 

```bash
juliaup add 1.10.5
juliaup default 1.10.5
```

for further details concerning `juliaup` go here:

[https://github.com/JuliaLang/juliaup?tab=readme-ov-file#using-juliaup]

### cloning the sga_index_consensus repository

Now that the dependencies are setup we clone the repository

```bash
cd ~
git clone https://github.com/Mullinslab/sga_index_consensus.git
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

### Workflow

Snakemake makes use of a `Snakefile` and a `config` file to specify 
input and output files for each rule and to set parameters for each 
rule in the workflow. Global parameters that can be changed by the 
user editing the `Snakefile` are as follows:

```
# sga-index-consensus parameters
# demux
chunk_size = 100000  		# default 100000
index_type = "Index_primer" # default "Index_Primer"
error_rate = 0.01    		# default 0.01
min_length = 2000    		# default 2100
max_length = 6000    		# default 4000
max_reads = 1000       		# default 1000 reads per sample,
                         	# use something large for no downsampling
# consensus
min_reads = 5				# default 5
# contam
cluster_thresh = 0.015   	# default 0.015
proportion_thresh = 0.2  	# default 0.2
dist_thresh = 0.015      	# default 0.015
# postproc
agreement_thresh = 0.7   	# default 0.7 for PacBio reads, 0.65 for Nanopore reads
max_alignment_reads = 1100  # default 1100
                            # be sure value is ~10% larger than max_reads to avoid conflicts
```

Note that with the advent of PacBio Revio sequencer, the number of reads
per sample has grown to outstrip memory available on standard CPUs. 
To enable a trouble free pipeline run, we now allow the user to specify
the maximum number of reads per sample using the `max_reads` parameter above.
Samples with reads exceeding this limit are then randomly sub-sampled to
reduce the number of reads accordingly.

If read counts are very high, alignments created for read collections
can take a long time, `max_alignment_reads` allows the user to create
alignements only for read collections below a specified count. 

### Configuration

To configure the workflow, use the format below to define
your config construction. 
It should follow the same format as shown below:

```yaml
Dataset1:
  Sample1:
    fwd_index: Index_F21
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    fwd_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample2:
    fwd_index: Index_F22
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    fwd_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample3:
    fwd_index: Index_F23
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    fwd_primer: GATTGTGTGGCARGTAGACAGRATG
```

The primer sequences provided will be used for demultiplexing and will be trimmed
from the final sequences. **fwd_primer** and **rev_primer** are the PCR primers 
used in the 2nd rd of PCR. **fwd_index** and **rev_index** are the index primers used to 
tag each sample. Index primer set is defined in the file Index_primer_sheet.csv.

To produce the config file required to run the pipeline for any given
sample set, a csv must first be created with the Index and PCR primer
information for each sample. It is critical to name this file correctly
as the name of this csv file is copied in as the dataset name in the
config file. This dataset name will need to match the CCS read file
(fastq.gz) and will be used to name the output data from the pipeline.

Example: Dataset1.csv goes with dataset Dataset1.fastq.gz

An example primer sheet can be found in the bash-scripts folder. (Dataset1.csv)

```
bash
cd bash-scripts/
python3 primersheet2config_SGA-Index.py Dataset1.csv
```

After creating the config file the contents must be copied into the
existing config.yaml file.

This config.yaml file is the main config file used per pipeline run,
dictating the dataset name, what samples are included, and also any
primer sequences and panel files required for each sample. This is
created from the primersheet2config_SGA-Index.py file, as described
above, and moved to the sga-index-consensus folder or copied to the
existing config.yaml file. Double check this config file looks correct
before a run to ensure primers are correct. A sample which gives an
issue or crashes the pipeline can be commented out to prevent the
pipeline from attempting to run it.

The snakefile contains the rules for the pipeline and also parameters
that may need to be changed depending on the run. For example, with
recent NFLG samples, the max_length parameter was changed to 10,000
(from the default 4000) and min_length changed to 10 (from the default
2100). This file mostly remains the same, other than any specific
snakefile parameter changes. Be sure that the snakefile is pointing to
the correct config file (line 3, `configfile: "config.yaml"`) or that
the correct contents have been copied to the config.yaml file.

Raw CCS .fastq files should be placed in the `raw-reads/` subdirectory and named 
according to the the dataset name used in the `config.yaml` file, ie, `Dataset1.fastq`
for the example dataset.

The contamination check will compare all sequences to one another, but if it is desired for
another sequence (for example a lab strain) to be included any sequence can be added to the
`contam_panel.fasta` file in the panel folder. The sequence must be trimmed to the same length
as the target sequence it is to be compared to. Note that the contamination check will flag
sequences that are within the distance threshold `dist_thresh` specified in the snakefile
(default 1.5%).

### Preview and Execution

Preview jobs with Snakemake and run with {n} cores.

```
bash
#preview jobs
snakemake -np

#run
snakemake -j{n}
```

For more info on Snakemake, see https://snakemake.readthedocs.io/en/stable/

### Pipeline Results
The output files from the pipeline are all written to dataset/[dataset
name]. There are 4 summary reports and 3 folders containing output from
the pipeline. The reports are described below:

reject_report.csv - the number of reads rejected by the filter for CCS
read quality 

quality_report.csv - a table summarizing total reads,
quality reads, and total number of reads assigned to samples

demux_report.csv- the number of reads matching the index and PCR primers
for each sample in the config and how many reads remained after downsampling.

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

