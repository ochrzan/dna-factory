# PopFactory (Population Factory)

PopFactory is a tool for generating simulated genetic population data for use in genetic analysis tools.  PopFactory
 creates [VCF files](https://samtools.github.io/hts-specs/) that are zipped in bgzf format for easy use with samtools
 .  In order to generate large VCF files, PopFactory can scale out to multiple processes and horizontally to
  multiple nodes. 
 
 ## Installing
 
 PopFactory is written in Python 3. You will need python 3 installed and also pip. The python package requirements are
  all in the requirements.txt file and can usually be installed by: 
 ```
pip3 install -r requirements.txt
```
 
 ## Generating a Population
 
 To generate a VCF file, run pop_factory.py from the command line. Here is a typical scenario:
<pre>
python3 pop_factory.py -s <i>num_cases</i> -c <i>num_controls</i> -x <i>num_snps</i> -f <i>min_minor_allele_freq</i> 
</pre>
For example, the following command will generate a VCF file with 200 samples (100 cases with mutations configured in
 the pathogens.yml file and 100 controls), 100000 SNPs, and all selected SNPs will have a minor allele frequency >= 1%.
```shell script
python3 pop_factory.py -s 100 -c 100 -x 100000 -f 0.01
```
There are options to control the output location (default is a timestamp subdir in the populations directory), number
 of worker processes, compression level, and odds of a population member being male (having a Y chromosome). For all
  available options run:
  ```shell script
python3 pop_factory.py -h
```
### Output Files
PopFactory outputs a collection of files each time it is run.
* population.vcf.gz - VCF file zipped with bgzip
* population.fam - [fam file](https://www.cog-genomics.org/plink/1.9/formats#fam) with information on each sample. It
 will have the following columns:
 1. Family ID ('FID') - all samples with have different family IDs
 2. Within-family ID ('IID') - same as sample id in VCF header
 3. Within-family ID of father - always zero
 4. Within-family ID of mother - always zero
 5. Sex code ('1' = male, '2' = female, '0' = unknown)
 6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
 * pathogens.json - JSON file with the selected pathogenic/causal mutations
 * 
 
 ## Generating Very Large Files
 If you want to generate a very large VCF file, there are two ways to speed things up:
 1. Increasing the number of workers
 2. Running multiple jobs in parallel and merging the resulting VCF files

 ### Increasing Workers (Processes)
 PopFactory uses a parallel processing approach where there is one write process and multiple worker processes that
  feed rows to the writer process. You can increase the throughput of PopFactory by increasing the number of workers:
<pre>
python3 pop_factory.py -s 10000 -c 10000 -x 5000000 -f 0.01 <strong>-n 7 -z 2</strong>
</pre>
 The *-n* option in this case is creating 7 worker processes (there will be 8 in total with the writer process) and
  is setting a gzip compression level of 2 (faster writes). This is a suitable configuration for an 8 core server. It
   is recommended to never set the number of workers higher than the number of available cores - 1.
### Running Multiple Parallel Jobs
If you desire to make a very large VCF file (millions of SNPs), it may be desirable to generate portions of the file
 in parallel and merge them together. PopFactory supports easy use of [bcftools](http://www.htslib.org/doc/bcftools.html) 
 to index and merge generated VCF files.
 
 To run multiple jobs in parallel, do the following steps.
 1. Run the tool once 
 ## Pathogens.yml (pathogens config)
 
 ## Downloading RefSNP Data
 
 ### Database Config
 
 

