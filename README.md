# StrainIQ - Strain *I*dentification and *Q*uantification
StrainIQ is an n-gram based method to identify and quantify microbial communities in metagenomics samples.

**Dependencies:**
* Python >= 3.6.2
* BioPython 1.72
* Pandas 0.23.4

**Usage**
```
StrainIQ.py [-h] -p Program [-n n-size]
                   [-glist file with a list of reference genome]
                   [-ng total number of genomes in the dasem] [-dsem DSEM]
                   [-sample sample] [-prediction -prediction]

optional arguments:
  -h, --help            show this help message and exit
  -p Program            Program to run. Example: builder, identifier,
                        quantifier
  -n n-size             Size of the n-gram
  -glist file with a list of reference genome
                        Tab delimited file with list of genomes to be included
                        in the DSEM. Format: gid filelocation
  -ng total number of genomes in the dasem
                        The number of genomes in the DSEM. This parameteris
                        needed if the DSEM is rebuilt using different sets for
                        references.
  -dsem DSEM            The DSEM for this bodysite.
  -sample sample        Input metagenomic sample. Convert the reverse reads in
                        case of paired end sequencing before bombining the
                        reads together
  -prediction -prediction
                        The prediction file produced by the identifier
```

StrainIQ has three main parts: Builder, Identifier, and Quantifier. 

**StrainIQ Builder:** StrainIQ-Builder generates a DSEM for a body site. This is run only once for each body site at the front end of the process and repeated as and when the genomes need to be updated in a DSEM. It takes n-size and the list of genomes in a body site as input. The –glist is a tab-delimited file with user-assigned genome-id and genome file location. Along with the output DSEM, the builder also creates a configuration file for use with identification and quantification steps. These files are available for [download](https://drive.google.com/drive/folders/18uutqK1ExeYYCQFdh8RM3PyMZS-pO4MA?usp=sharing) for GI Tract, Blood and Urogenital tract.

**DSEM Building**

```
python StrainIQ.py 
	    -p builder # Use builder program to build DSEM
	    –n 21 # Default n-gram size for GI.
	    –glist <genome list> # genome file and location.
```

**StrainIQ Identifier:** StrainIQ-Identifier takes –dsem and sample name as input to identify the microbes in the given sample. The identifier refers to the configuration file for the additional parameters. The configuration file is generated as a part of the DSEM building and has the same name as DSEM with the .conf extension. I addition to DSEM and configuration file, the identifier also refers to the Map file and Taxonomy file.

```
python StrainIQ.py 
    -p identifier # Use identifier program for identification
    -dsem gi.dsem # Choose appropriate DSEM for the body site
    -sample sample1.fastq #Provide sample in fastq format
```
**StrainIQ Quantifier:** StrainIQ-Quantifier takes –dsem, -sample and –prediction file as input to calculate the abundance of microbes in the metagenomic sample. 
```
python StrainIQ.py 
    -p quantifier #Use quantifier program for quantification
    -dsem gi.dsem #Choose appropriate DSEM for the body site
    -sample sample1.fastq # sample in fastq format
    -prediction sample1.prediction # list identified genomes
```
