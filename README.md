# PopFactory (Population Factory)

PopFactory is a tool for generating simulated genetic population data for use in genetic analysis tools (especially
 GWAS tools).  PopFactory creates [VCF files](https://samtools.github.io/hts-specs/) that are zipped in bgzf format
  for easy use with samtools.  In order to generate large VCF files, PopFactory can scale out to multiple processes and horizontally to
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
 the deleterious.yml file and 100 controls), 100000 SNPs, and all SNPs in the file will have a minor allele frequency
  greater than 1%.
```shell script
python3 pop_factory.py -s 100 -c 100 -x 100000 -f 0.01
```
There are options to control the output location (default is a new subdir with a timestamp name in the populations
 direcotry), number
 of worker processes, compression level, and odds of a population member being male (having a Y chromosome). For all
  available options run:
  ```shell script
python3 pop_factory.py -h
```
### Output Files
PopFactory outputs a collection of files each time it is run.
* population.vcf.gz - VCF output file zipped with bgzip. Controls will follow the regular SNP frequency distribution
 for SNPs in the file.  Cases (population memeber with an affliction) will follow regular SNP frequencies except for
  selected deleterious SNPs.
* population.fam - [fam file](https://www.cog-genomics.org/plink/1.9/formats#fam) with information on each sample. It
 will have the following columns:
 1. Family ID ('FID') - all samples with have different family IDs
 2. Within-family ID ('IID') - same as sample id in VCF header
 3. Within-family ID of father - always zero
 4. Within-family ID of mother - always zero
 5. Sex code ('1' = male, '2' = female, '0' = unknown)
 6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
 * deleterious.json - JSON file with the selected deleterious SNPs. This file can be used as an input
   to future PopFactory runs.
 * pop_deleterious.txt - text file which maps sample IDs to deleterious groups and deleterious mutations
 * snps.json.gz - JSON file with the SNPs that were used to generate the VCF file. This file can be used as an input
  to future PopFactory runs.
 
 ## Generating Very Large Files
 If you want to generate a very large VCF file, there are two ways to speed things up:
 1. Increase the number of workers
 2. Run multiple jobs in parallel and merge the resulting VCF files

 ### Increasing Workers (Processes)
 PopFactory uses a parallel processing approach where there is one write process and multiple worker processes that
  feed rows to the writer process. You can increase the throughput of PopFactory by increasing the number of workers:
<pre>
python3 pop_factory.py -s 10000 -c 10000 -x 5000000 -f 0.01 <strong>-n 7 -z 2</strong>
</pre>
 The *-n* option in this case is creating 7 worker processes (there will be 8 processes in total when you include the
  writer process). This command also uses the **-z** option to set a gzip compression level of 2 (faster writes
  /less compression). This is a suitable configuration for an 8 core server. It
   is recommended to never set the number of workers higher than the number of available cores - 1.
### Running Multiple Parallel Jobs
If you desire to make a very large VCF file (millions of SNPs), it may be desirable to generate portions of the file
 in parallel and merge them together. PopFactory supports easy use of [bcftools](http://www.htslib.org/doc/bcftools.html) 
 to index and merge generated VCF files.
 
 To run multiple jobs in parallel, perform the following steps.
 <ol><li> Run the tool once to generate/select SNPs and deleterious SNPs
<pre>
python3 pop_factory.py -s 0 -c 0 -p "path_to_deleterious.yml" -f 0.01 -x 10000000 --outdir "path_to_output_dir"
</pre>
This will select 10 million SNPs and pick deleterious SNPs (based on the yml config), but will make an empty VCF file.
</li>
<li>Run the tool again and use the outputs from step 1 as inputs to generate a population
<pre>
python3 pop_factory.py -s 1000 -c 1000 --outdir "path_to_output_dir2" <b>--snps_file
 "path_to_output_dir"/snps.json.gz --deleterious_file "path_to_output_dir"/deleterious.json</b>
</pre>
This will make a VCF file with 1000 cases and 1000 controls using the passed in SNPs and deleterious SNPs
</li>
<li>Run the command as many times as you want using incremental offsets to make sure all samples get unique IDs. The
 offset number is added to the generated sample IDs. You will want to use a different offset for each file to ensure
  the sample ID space is unique.
<pre>
python3 pop_factory.py  <b>--offset 1000</b> -s 1000 -c 1000 --outdir "path_to_output_dir3" --snps_file
 "path_to_output_dir"/snps.json.gz --deleterious_file "path_to_output_dir"/deleterious.json
</pre>
</li>
<li>Use bcftools to merge files as desired.
<pre>
bcftools index "path1/population.vcf.gz"
bcftools index "path2/population.vcf.gz"
bcftools merge "path1/population.vcf.gz" "path2/population.vcf.gz"
</pre>
</li>
</ol>

## Deleterious.yml (deleterious config)

PopFactory selects deleterious mutations based on a configuration file in yaml format. The default location is
 "deleterious.yml" in the install directory. You can modify this file or create your own deleterious.yml config files
  to mimic many hypothetical scenarios. You can create multiple deleterious groups, polygenic groups, rare or common
   groups
 , etc. 
  
  Here are the options:
 
 ```shell script
group_name:
  mutation_weights:
    - 1
    - 0.5
    - 0.2
  num_instances: 1
  population_weight: 2
  max_minor_allele_freq: 0.10
  min_minor_allele_freq: 0.01
```
* group_name - Can be any string to name the deleterious group
* mutation_weights - list of items with the weights of each deleterious mutation. The number of items in list
 determines the number of SNPs selected and the weights determine how many SNPs each afflicated case has. Every case
  will have mutations (minor alleles) that add up to a weight >= 1. 
  
  Example 1: 4 SNPs, each case has 2 of 4
 ```yaml
mutation_weights:
  - 0.5
  - 0.5
  - 0.5
  - 0.5
```  
Example 2: 5 SNPs, the first SNP has double the weight. 3 SNPs needed unless the 1st SNP is selected for the case.
```yaml
mutation_weights:
  - 0.7
  - 0.35
  - 0.35
  - 0.35
  - 0.35
```
Example 3: 3 SNPs, 1 out of 3 is needed
```yaml
mutation_weights:
  - 1
  - 1
  - 1
```
* num_instances - number of copies of this group (SNPs are reselected)
* population_weight - Proportional weight this group has in the overall afflicted population. For example, if the sum
 of all population_weights is 20 and this group has a weight of 2, it will exist in 10% of cases.
* max_minor_allele_freq - highest minor allele frequency allowed for SNPs in this group
* min_minor_allele_freq - lower bound for minor allele frequency for SNPs in this group
## Using RefSNP Data

By default, PopFactory generates SNPs based on a continuous distribution function created from RefSNP frequency data
 (see snp_freq_cdf.csv for values). You can also download real RefSNP data from NIH and use it for SNP selection. The
  RefSNP dataset is very large and building a RefSNP database can take some time. The downloaded data is minimized to
  only include data useful for PopFactory and stored in a SQL database (default SQLite). A fullly downloaded RefSNP
   SQLite DB can
  be over 70 GB on disk.
  
  Here is how to download the data and build a DB:
  
```shell script
python3 download.py
```
If you want to store the data in a different database than the default, supply a [sqlalchemy connection string](https://docs.sqlalchemy.org/en/13/core/engines.html) in **db
.yml**.

You can also download partial data (or fill in missing data) on a chromosome by chromosome basis.
```shell script
python3 download.py -c 1,2,X -n 4
```
This will download data for chromosomes 1, 2, and X. It will use 4 worker processes for downloading. It will only
 delete data for the specified chromosomes and will leave the other chromosome data in place. This can be useful if
  the download job fails or needs to be stopped. The download process also used MD5 hashing to avoid redownloading
   files from NIH which have not changed (assuming the local tmp_download files are still on disk).

To then use the built RefSNP DB, use the **-l** flag when running PopFactory. This will cause the code to pull SNP
 data from the DB instead of generating simulated SNPs. For example:
 ```shell script
python3 pop_factory.py -l -s 100 -c 100 -x 10000 -f 0.02
``` 
 

