"""
Generates fake data for similating possible scenarios for use in PLINK.

1. Read config/command
2. Read SNP data
3. Generate PED/MAP files based on command and SNP data

Similar to http://cnsgenomics.com/software/gcta/#GWASSimulation ?
"""

import getopt
import sys
import glob
import yaml
from common.snp import RefSNP, Allele, is_haploid
from download import OUTPUT_DIR
import re
import numpy
import os
from datetime import datetime
import gzip

MIN_SNP_FREQ = 0.005
OUTPUT_DIR = "populations"
SNP_DIR = "output"
TEST_MODE = False

class SNPTuples:
    """ Class for holding compressed snp and probability data
    """

    def __init__(self, snp_id):
        self.id = snp_id
        self.tuples = []

    def add_tuple(self, inserted, range_end):
        self.tuples.append((inserted, range_end))

    def pick_snp_value(self, random_roll):
        for nt_letter, prob in self.tuples:
            if prob > random_roll:
                return nt_letter

    def pick_pathogen_value(self):
        """
        Returns the second most probable snp value as the "pathogen" mutation
        :return: second most probable nucleotide for this snp
        """
        if len(self.tuples) == 1:
            return self.tuples[0][0]
        prev_prob = 0
        sorted_probabilities = []
        for nt_letter, prob in self.tuples:
            sorted_probabilities.append((nt_letter, prob - prev_prob))
            prev_prob = prob
        sorted_probabilities.sort(key=lambda x: x[1], reverse=True)
        return sorted_probabilities[1][0]

class PopulationFactory:

    # number of subgroups with phenotype, total number of hidden mutations
    def __init__(self, num_groups, num_mutations):
        self.num_groups = num_groups
        self.num_mutations = num_mutations
        self.pathogens = {}
        self.ordered_snps = []
        self.snp_count = 0
        self.population_dir = OUTPUT_DIR

    def generate_population(self, control_size, test_size, male_odds, snp_data):
        """Generate a simulated population based on the number of groups, mutations, size of test group,
        size of control group and the snp dictionary.
        1. Determine which snps will be the hidden pathogenic snps
        2. Generate control data using random generated population
        3. Generate test data based on hidden pathogens and random otherwise
        Use numpy to generate random number en mass
        """
        subdir = datetime.now().strftime("%Y%m%d%H%M")
        numpy.random.seed(int(datetime.now().strftime("%H%M%S")))
        self.population_dir = OUTPUT_DIR + "/" + subdir + "/"
        os.makedirs(self.population_dir, exist_ok=True)

        self.pick_pathogen_snps(snp_data)
        # Create map file

        self.output_map_file(snp_data)
        # Create control population
        self.output_population(control_size, True, male_odds)
        self.output_population(test_size, False, male_odds)
        return

    def output_map_file(self, snp_data):
        """
        Outputs a map file used by plink (snps). Also sets self.ordered_snps to be a list of lists of tuples.
        One list per SNP that has a tuple per allele with the inserted value and probability range. For instance
        if a SNP has 3 alleles A (55%), T (25%), C (20%) the tuples would be ("A",0.55), ("T",0.8), ("C", 1.0)
        :param snp_data: catalog of ReFSNP data (as loaded from yml files)
        :return:
        """

        with open(self.population_dir + "population.map", 'w') as f:
            for chromo, snps in snp_data.items():
                # TODO save Y chromosome and MT DNA separtely.
                chromo_snps = []
                for snp in snps.values():

                    f.write("%s\trs%s\t0\t%s\n" % (chromo, snp.id, list(snp.alleles.values())[0].position))
                    # Make numpy array in same order for reference
                    # ideally each item would have a tuple of inserted value and probability
                    running_allele_count = 0
                    snp_tuple = SNPTuples(snp.id)
                    for allele in snp.alleles.values():
                        snp_tuple.add_tuple(allele.inserted,
                                           (allele.allele_count + running_allele_count) / allele.total_count)
                        running_allele_count += allele.allele_count
                    chromo_snps.append(snp_tuple)
                # Save each chromosome separately, but in an ordered list of tuples so the line up with the map file
                self.snp_count += len(chromo_snps)
                self.ordered_snps.append((chromo, chromo_snps))

    def output_population(self, size, is_control, male_odds):
        """
        Output a population file of two nucleotide values per SNP. Correctly outputs duplicate values if the
        chromosome is haploid
        :param size: size of the generated population
        :param is_control: control population with no hidden pathogens
        :param male_odds: odds of a person being a biological male
        :return:
        """
        with open(self.population_dir + "population.ped", 'a+') as f:
            for i in range(size):
                # Roll the dice for each snp and each allele. This will be a bit long for boys, but will work
                randoms = numpy.random.rand(self.snp_count * 2)
                is_male = randoms[0] < male_odds
                snp_values = []
                i = 0
                row = 1
                for chromo, snps in self.ordered_snps:
                    if not is_male and chromo == 'Y':
                        continue  # Skip Y snps for women
                    for snp in snps:
                        if is_control or not self.is_pathogen(snp.id):
                            random_roll = randoms[i]
                            i += 1
                            selected_nt = snp.pick_snp_value(random_roll)
                            if is_haploid(chromo, is_male):
                                other_nt = selected_nt
                            else:
                                random_roll = randoms[i]
                                i += 1
                                other_nt = snp.pick_snp_value(random_roll)
                            snp_values.append(selected_nt)
                            snp_values.append(other_nt)
                        else:
                            selected_nt = snp.pick_pathogen_value()
                            snp_values.append(selected_nt)
                            snp_values.append(selected_nt)
                            # TODO make it so pathogens can be recessive or dominant
                #Output row - Family ID, Indiv ID, Dad ID, MomID, Sex, affection, snps
                if is_male:
                    sex = 1
                else:
                    sex = 2
                if is_control:
                    affection = 1
                else:
                    affection = 2
                f.write("\t".join(map(lambda x: str(x), [row, row, 0, 0, sex, affection])) + "\t" + "\t".join(snp_values) + "\n")
                row += 1

    def pick_pathogen_snps(self, snp_data):
        """
        Pick and store the snps that are the pathogens. Randomly? pick num_mutations from the snps
        :param snp_data: SNP candidates
        :return:
        """
        # Kinda inefficient... Could just pick based on total size and use indexes, but that is more complex
        snp_id_list = []
        for chromo, snps in snp_data.items():
            snp_id_list.extend(list(map(lambda x: x.id, snps.values())))
        for snp_id in numpy.random.choice(snp_id_list,  self.num_mutations):
            self.pathogens[snp_id] = snp_id
        with open(self.population_dir + "pathogens.txt", 'w') as f:
            for pathogen_id in self.pathogens:
                f.write(str(pathogen_id) + "\n")

    def is_pathogen(self, snp_id):
        return snp_id in self.pathogens


def load_snps(dir, min_freq):
    # Seems the entire refSNP db might be in the order of 400 million SNPs so filtering will be needed
    # 95% of mutations in any persons's genome are from common mutations (>1% odds), though.
    # We may need to shrink to only include a subset using a min frequency threshold
    if TEST_MODE:
        snp_file_list = glob.glob(dir + "/*chr*sample.yml")
    else:
        snp_file_list = glob.glob(dir + "/*chr*.yml*")
    loaded_snps = {}
    for snp_file in snp_file_list:
        chrom_snps = {}
        chr_search = re.search('chr([0-9XYMT]+)', snp_file, re.IGNORECASE)
        if chr_search:
            chromosome = chr_search.group(1)
        else:
            chromosome = 'unknown'
        loaded_snps[chromosome] = chrom_snps
        open_fn = open
        if snp_file.endswith(".gz")
            open_fn = gzip.open
        with open_fn(snp_file) as f:
            items = yaml.load(f, Loader=yaml.FullLoader)
            for name, alleles in items.items():
                if not alleles:
                    continue
                # find the most common allele
                max_allele_count = 0
                total_count = 0
                is_valid_for_plink = True
                for allele in alleles:
                    if not allele['deleted'] or not allele['inserted'] or allele['total_count'] < 1000:
                        #Skip inserts, deletes and small samples
                        is_valid_for_plink = False
                        break
                    if len(allele['inserted']) > 1 or len(allele['deleted']) > 1:
                        #Skip multi-NT snps
                        is_valid_for_plink = False
                        break
                    if allele['allele_count'] > max_allele_count:
                        max_allele_count = allele['allele_count']
                        total_count = allele['total_count']
                if not is_valid_for_plink:
                    continue
                common_allele_freq = max_allele_count / total_count
                if common_allele_freq <= (1 - min_freq):
                    # If passes freq filter, then save it
                    snp = RefSNP(name)
                    for allele_attr in alleles:
                        allele = Allele(allele_attr['deleted'], allele_attr['inserted'],
                                        allele_attr['seq_id'], allele_attr['position'])
                        allele.allele_count = allele_attr['allele_count']
                        allele.total_count = allele_attr['total_count']
                        snp.put_allele(allele)
                    chrom_snps[snp.id] = snp
    return loaded_snps


def print_help():
    print("""
    Accepted Inputs are:
    -g number of hidden subgroups all with same phenotype
    -n total number of different mutations that can cause a subphenotype (min 1 per subgroup)
    -s size of test group (afflicted group)
    -c size of control group
    -f min frequency for a SNP to be included in the list of SNPs, default is 0.005
    """)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h?ngfs:c:t", ["help"])
    except getopt.GetoptError as err:
        print(err.msg)
        print_help()
        sys.exit(2)
    subgroups = 1
    num_mutations = 1
    min_freq = MIN_SNP_FREQ
    male_odds = 0.5
    for opt, arg in opts:
        if opt in ('-h', "-?", "--help"):
            print_help()
            sys.exit()
        elif opt in "-n":
            num_mutations = int(arg)
        elif opt in "-g":
            subgroups = int(arg)
        elif opt in "-s":
            size = int(arg)
        elif opt in "-c":
            control_size = int(arg)
        elif opt in "-f":
            min_freq = float(arg)
        elif opt in "-m":
            male_odds = float(opt)
        elif opt in "-t":
            TEST_MODE = True
    snps = load_snps(SNP_DIR, min_freq)
    pop_factory = PopulationFactory(subgroups, num_mutations)
    pop_factory.generate_population(control_size, size, male_odds, snps)

if __name__ == '__main__':
    main(sys.argv[1:])

