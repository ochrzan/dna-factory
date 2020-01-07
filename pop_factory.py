"""
Generates fake data for similating possible scenarios for use in PLINK.

1. Read config/command
2. Read SNP data
3. Generate PED/MAP files based on command and SNP data

"""

import getopt
import sys
import glob
import yaml
from common.snp import RefSNP, Allele
from download import OUTPUT_DIR
import re


MIN_SNP_FREQ = 0.005


class PopulationFactory:

    # number of subgroups with phenotype, total number of hidden mutations
    def __init__(self, num_groups, num_mutations):
        self.num_groups = num_groups
        self.num_mutations = num_mutations


    def generate_population(self, control_size, test_size, phenotype, snp_data):
        """Generate a simulated population based on the number of groups, mutations, size of test group,
        size of control group and the snp dictionary.
        1. Determine which snps will be the hidden pathogenic snps
        2. Generate control data using random generated population
        3. Generate test data based on hidden pathogens and random otherwise
        Use numpy to generate random number en mass
        Use numpy to generate snps output en mass?"""

        return


def load_snps(dir, min_freq):
    # TODO Seems the entire refSNP db might be in the order of 400 million SNPs.
    # 95% of mutations in any persons's genome are from common mutations (>1% odds), though.
    # We may need to shrink to only include a subset. Add a min frequency threshold
    snp_file_list = glob.glob(dir + "/*.yml")
    loaded_snps = {}
    for snp_file in snp_file_list:
        chrom_snps = []
        chr_search = re.search('chr([0-9XYMT]+)', snp_file, re.IGNORECASE)
        if chr_search:
            chromsome = chr_search.group(1)
        else:
            chromosome = 'unknown'
        loaded_snps[chromosome] = chrom_snps
        with open(snp_file) as f:
            items = yaml.load(f)
            for name, alleles in items.items():
                # TODO load in yml files, filter based on min freq, return snps
                # Assumes first allele is the more common allele
                common_allele_freq = alleles[0]['allele_count'] / alleles[0]['total_count']
                if alleles[0]['allele_count'] / alleles[0]['total_count'] < 0.5:
                    raise Exception('First Allele in SNP expected to have frequency greater than 50%.'
                                    ' Frequency is {} ({} / {})'.format(common_allele_freq,
                                                                        alleles[0]['allele_count'],
                                                                        alleles[0]['total_count']))
                if common_allele_freq <= (1 - min_freq):
                    snp = RefSNP(name)
                    for allele_attr in alleles:
                        allele = Allele(allele_attr['deleted'], allele_attr['inserted'],
                                        allele_attr['seq_id'], allele_attr['position'])
                        allele.allele_count = allele_attr['allele_count']
                        allele.total_count = allele_attr['total_count']
                        snp.put_allele(allele)
                    chrom_snps.append(snp)
    return loaded_snps




def print_help():
    print("""
    Accepted Inputs are:
    -s number of hidden subgroups all with same phenotype
    -n total number of different mutations that can cause a subphenotype (min 1 per subgroup)
     
    """)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h?ngfs:c:", ["help"])
    except getopt.GetoptError as err:
        print(err.msg)
        print_help()
        sys.exit(2)
    subgroups = 1
    num_mutations = 1
    min_freq = MIN_SNP_FREQ
    disease_name = 'HAS_DISEASE'
    for opt, arg in opts:
        if opt in ('-h', "-?", "--help"):
            print_help()
            sys.exit()
        elif opt in "-n":
            num_mutations = int(arg)
        elif opt in "-g":
            subgroups = int(arg)
        elif opt in "s":
            size = int(arg)
        elif opt in "c":
            control_size = int(arg)
        elif opt in "f":
            min_freq = float(arg)
    snps = load_snps(OUTPUT_DIR, min_freq)
    pop_factory = PopulationFactory(subgroups, num_mutations)
    results = pop_factory.generate_population(control_size, size, disease_name, snps)

if __name__ == '__main__':
    main(sys.argv[1:])

