"""
Generates VCF data for simulating possible scenarios for use in PLINK or other bio analysis tools.
"""
import gc
import json
import random
import sys
import argparse
import time
from multiprocessing import Process, Queue
import queue
import heapq
from Bio import bgzf
import numpy
import os
from datetime import datetime
import gzip
from yaml import load

from common.db import db
from common.timer import Timer
from common.snp import RefSNP, Allele, is_haploid, \
    CHROMOSOME_PROB, CHROMOSOME_LIST, CHROMOSOME_MAX_POSITION, stripe_list
from definitions import ROOT_DIR

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

MIN_SNP_FREQ = 0.005
MIN_TOTAL_COUNT = 1000
OUTPUT_DIR = os.path.join(ROOT_DIR, "populations")


def gen_vcf_header(fam_data):
    header = "##fileformat=VCFv4.3\n"
    header += "##filedate=%s\n" % datetime.now().strftime("%Y%m%d %H:%M")
    header += "##source=SNP_Simulator\n"
    header += '##FILTER=<ID=q10,Description="Quality below 10">\n'
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    header += "\t".join(map(lambda x: str(x.person_id), fam_data)) + "\n"
    return header


class SampleInfo:
    """
    Class for individual sample metadata for use in a .fam file. Also holds deleterious data for this individual
    """

    def __init__(self, family_id, person_id, father_id, mother_id, sex: int, is_control: bool, deleterious_snps: dict):
        assert person_id
        self.person_id = person_id
        self.family_id = family_id
        self.father_id = father_id
        self.mother_id = mother_id
        self.sex = sex
        self.is_control = is_control
        self.deleterious_snps = deleterious_snps

    def to_fam_format(self):
        if self.is_control:
            pheno_code = 1
        else:
            pheno_code = 2
        return "%i\t%i\t%i\t%i\t%i\t%i\t\n" % \
               (self.family_id, self.person_id, self.father_id, self.mother_id, self.sex, pheno_code)

    def is_male(self):
        return self.sex == 1


class SNPTuples:
    """ Class for holding compressed snp and probability data
    """

    def __init__(self, snp_id, chromosome, position):
        self.id = snp_id
        self.chromosome = chromosome
        self.tuples = []
        self.position = position

    def add_tuple(self, inserted, range_end):
        self.tuples.append((inserted, range_end))

    def pick_snp_value(self, random_roll):
        for nt_letter, prob in self.tuples:
            if prob > random_roll:
                return nt_letter

    def pick_allele_index(self, random_roll):
        for i, tupl in enumerate(self.tuples):
            if tupl[1] >= random_roll:
                return i

    def minor_allele_tuple(self):
        """
        Returns the second most probable allele and it's frequency
        :return: second most probable nucleotide for this snp and it's frequency
        """
        return self.tuples[1]

    def ref_allele_tuple(self):
        """
        Returns the most probable allele and it's frequency
        :return: second most probable nucleotide for this snp and it's frequency
        """
        return self.tuples[0]

    def alt_alleles(self):
        if len(self.tuples) == 1:
            return self.tuples[0][0]
        if len(self.tuples) == 2:
            return self.tuples[1][0]
        return ",".join(map(lambda x: x[0], self.tuples[1:]))

    def __str__(self):
        json_hash = {"id": self.id, "chromosome": self.chromosome, "position": self.position}
        if len(self.tuples) > 0:
            json_hash["tuples"] = {}
        for t in self.tuples:
            json_hash["tuples"][t[0]] = t[1]
        return json.dumps(json_hash)

    @classmethod
    def from_json(cls, json_line):
        ref_obj = json.loads(json_line)
        snp_tuples = cls(ref_obj['id'], ref_obj["chromosome"], ref_obj["position"])
        if "tuples" in ref_obj:
            for i, f in ref_obj['tuples'].items():
                snp_tuples.add_tuple(i, f)
        return snp_tuples


class SnpFactory:

    def __init__(self, cdf_matrix):
        """
        Expects a numpy array with col 1 the maf and col2 the CDF. Easiest to use init_from_cdf_file to load
        values from a csv.
        """
        self.sorted_maf = cdf_matrix[:, 0]
        self.cdf = cdf_matrix[:, 1]
        cdf_shift = numpy.roll(numpy.append(self.cdf, 0), 1)
        self.pdf = self.cdf - cdf_shift[0:len(self.cdf)]

    @classmethod
    def init_from_cdf_file(cls, file="snp_freq_cdf.csv"):
        """
        Inits from a csv file that maps snp MAF frequency to it's CDF (cumulative distribution function).
        Default is to load from a CSV created from RefSNP freq probabilities across the entire genome.
        :param file: csv file that maps a snp MAF freq to it's CDF. Col 1 is the MAF, Col 2 is CDF function.
        Expects a header row
        :return: new SnpFactory object
        """
        freq_counts = numpy.loadtxt(os.path.join(ROOT_DIR, file), skiprows=1, delimiter=",")
        return cls(freq_counts)

    def gen_mafs(self, size, min_maf):
        start_maf_index = 0
        for i, m in enumerate(self.sorted_maf):
            if min_maf <= m:
                start_maf_index = i
                break
        return numpy.random.choice(self.sorted_maf[start_maf_index:], size=size,
                                   p=self.pdf[start_maf_index:] * 1 / numpy.sum(self.pdf[start_maf_index:]))

    def gen_chromosomes(self, size):
        return numpy.random.choice(CHROMOSOME_LIST, size=size, p=CHROMOSOME_PROB)

    def random_snp_tuples(self, size, min_maf=0.005):
        """
        Generate a random SNPTuples object with a MAF based on picking a random chromosome, position and an
        MAF that fits the RefSNP db distribution function for reported SNPs. Min MAF is 0.005. Max MAF is 0.50
        MAFs are stepwise with 0.005 increments.
        :param size: number of snptuples to generate
        :return: Generated SNPTuples
        """
        chromosomes = self.gen_chromosomes(size)
        mafs = self.gen_mafs(size, min_maf)
        position_randoms = numpy.random.random(size)
        nt_randoms = numpy.random.choice(["A", "T", "C", "G"], size=size)
        snp_tuples_list = []
        for n, c in enumerate(chromosomes):
            snp_tuple = SNPTuples(n + 1, c, int(position_randoms[n] * CHROMOSOME_MAX_POSITION[c]))
            snp_tuple.add_tuple(nt_randoms[n], 1 - mafs[n])
            remaining_nt = ["A", "T", "C", "G"]
            remaining_nt.remove(nt_randoms[n])
            alt_allele = random.choice(remaining_nt)
            snp_tuple.add_tuple(alt_allele, 1.0)
            snp_tuples_list.append(snp_tuple)
        return snp_tuples_list


class PopulationFactory:

    # number of subgroups with phenotype, total number of hidden mutations
    def __init__(self, num_processes=1, generate_snps=False, male_odds=0.5, deleterious_config=None,
                 deleterious_list_path=None, sample_id_offset=0, snps_path=None, output_path=None):
        self.deleterious = {}
        self.ordered_snps = []
        self.snp_count = 0
        if output_path:
            self.population_dir = output_path
            if not self.population_dir.endswith(os.path.sep):
                self.population_dir += os.path.sep
        else:
            subdir = datetime.now().strftime("%Y%m%d%H%M")
            self.population_dir = os.path.join(OUTPUT_DIR, subdir) + os.path.sep
        self.male_odds = male_odds
        if num_processes > 0:
            self.num_processes = num_processes
        else:
            self.num_processes = 1
        self.generate_snps = generate_snps
        self.deleterious_config = deleterious_config
        if sample_id_offset:
            self.sample_id_offset = sample_id_offset
        else:
            self.sample_id_offset = 0
        self.deleterious_list_path = deleterious_list_path
        self.snps_path = snps_path

    @Timer(logger=print, text="Finished Generating Population in {:0.4f} secs.")
    def generate_population(self, control_size, test_size, min_freq, max_snps,
                            compression_level=6):
        """Generate a simulated population based on the number of groups, mutations, size of test group,
        size of control group and the snp dictionary.
        1. Determine which snps will be the hidden deleterious snps
        2. Generate control data using random generated population
        3. Generate test data based on hidden deleterious and random otherwise
        Uses numpy to generate random numbers en mass
        """
        numpy.random.seed(int(datetime.now().strftime("%H%M%S")))
        os.makedirs(self.population_dir, exist_ok=True)
        if self.snps_path:
            self.load_snps_file()
        else:
            if self.generate_snps:
                snp_factory = SnpFactory.init_from_cdf_file()
                self.ordered_snps = snp_factory.random_snp_tuples(max_snps)
            else:
                self.load_snps_db(min_freq, max_snps)
        self.ordered_snps.sort(key=lambda x: (x.chromosome, x.position))
        if not self.snps_path:
            self.output_snps()
        gc.collect()
        if self.deleterious_list_path:
            self.load_deleterious()
        else:
            self.pick_deleterious_snps(self.ordered_snps, self.deleterious_config)

        # Create control population
        self.output_vcf_population(control_size, test_size, self.male_odds, compression_level)
        return

    @Timer(logger=print, name="PopFactory.output_snps()", text="Time to write snps file {:0.4f} seconds")
    def output_snps(self):
        with gzip.open(self.population_dir + "snps.json.gz", 'wt', compresslevel=5) as f:
            for t in self.ordered_snps:
                f.write(str(t) + "\n")

    def load_snps_file(self):
        """
        Load SNPTuples from a saved json file.
        :return:
        """
        with gzip.open(self.snps_path, 'rt') as f:
            for line in f:
                self.snp_count += 1
                self.ordered_snps.append(SNPTuples.from_json(line))

    def load_snps_db(self, min_freq, max_snps):
        """
        Load snps from DB and store as SNPTuples. Also output map file for plink.
        :param max_snps: Max number of snps to load
        :param min_freq: min Minor Allele frequency
        :return:
        """

        invalid_count = 0
        snps_result = db.connection.execute(
            "Select r.id, chromosome, maf, total_count,  deleted, inserted, position, allele_count "
            "from ref_snps r  "
            "join alleles a on r.id = a.ref_snp_id "
            "and r.maf >= %f and r.total_count >= %i" % (min_freq, MIN_TOTAL_COUNT)
        )
        current_snp_id = -1
        snp = None
        for snp_row in snps_result:
            if snp_row["id"] != current_snp_id:
                if snp and snp.valid_for_plink():
                    if self.snp_count >= max_snps - 1:
                        print("Hit max_snps size of %i. Stopping loading snps." % max_snps, flush=True)
                        break
                    self.add_snp_tuple(snp)
                    if self.snp_count % 100000 == 0:
                        print("Loaded %i snps. %s" % (self.snp_count, datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                              flush=True)
                else:
                    invalid_count += 1
                # otherwise new snp row
                snp = RefSNP.from_row_proxy(snp_row)

            # Added joined allele data every time
            snp.put_allele(Allele.from_row_proxy(snp_row))
            current_snp_id = snp_row["id"]
        self.add_snp_tuple(snp)
        print("Skipped Invalid:        %i" % invalid_count, flush=True)
        print("Total Loaded:           %i" % len(self.ordered_snps), flush=True)

    def add_snp_tuple(self, snp):
        """
        Adds snp to self.ordered_snps.
        self.ordered_snps is in the same order as the map file.
        One list per SNP that has a tuple per allele with the inserted value and cummulative distribution range.
        For instance if a SNP has 3 alleles A (55%), T (25%), C (20%) the tuples would be
         ("A",0.55), ("T",0.8), ("C", 1.0)"""
        running_allele_count = 0
        snp_tuple = SNPTuples(snp.id, snp.chromosome, snp.alleles[0].position)
        snp.alleles.sort(key=lambda x: x.allele_count, reverse=True)
        # Insert tuples in sorted order by frequency desc
        for allele in snp.alleles:
            snp_tuple.add_tuple(allele.inserted,
                                (allele.allele_count + running_allele_count) / snp.total_count)
            running_allele_count += allele.allele_count

            # Save each chromosome separately, but in an ordered list of tuples so the line up with the map file
        self.snp_count += 1
        self.ordered_snps.append(snp_tuple)

    @classmethod
    def pick_deleterious_groups(cls, deleterious_groups, pop_size):
        return random.choices(
            population=list(deleterious_groups),
            weights=list(map(lambda x: x.population_weight, deleterious_groups)),
            k=pop_size
        )

    def generate_fam_file(self, control_size, test_size, male_odds, deleterious_group_list):
        """

        :param control_size:
        :param test_size:
        :param male_odds:
        :param deleterious_group_list:
        :return: Data for each sample
        """
        control_id = 100000 + self.sample_id_offset
        test_id = 500000 + self.sample_id_offset
        randoms = numpy.random.rand(control_size + test_size)
        sample_data = []
        with open(self.population_dir + "population.fam", 'w') as f, \
                open(self.population_dir + "pop_deleterious.txt", "w") as pp:
            j = 0
            for i in range(control_size + test_size):
                is_control = i < control_size
                if is_control:
                    control_id += 1
                    iid = control_id
                else:
                    test_id += 1
                    iid = test_id
                if randoms[i] <= male_odds:
                    sex_code = 1
                else:
                    sex_code = 2

                if not is_control:
                    deleterious_group = deleterious_group_list[j]
                    j += 1
                    deleterious_snps = deleterious_group.select_mutations()
                    pp.write("%i\t%s\t" % (test_id, deleterious_group.name) +
                             "\t".join(map(lambda x: "rs" + str(x), deleterious_snps.keys())) + "\n")
                else:
                    deleterious_snps = None
                sample = SampleInfo(i + 1 + self.sample_id_offset * 2,
                                    iid, 0, 0, sex_code, is_control, deleterious_snps)
                sample_data.append(sample)
                f.write(sample.to_fam_format())

        return sample_data

    def output_vcf_population(self, control_size, test_size, male_odds, compression_level):
        """
        Output a population .vcf file and companion .fam file.
        :param test_size: size of control group
        :param control_size: size of cases/test group
        :param male_odds: odds of a person being a biological male
        :return:
        """

        if not self.ordered_snps:
            raise Exception("No SNPs to Process! Exiting.")
        # pick deleterious groups for population size
        deleterious_group_list = PopulationFactory.pick_deleterious_groups(list(self.deleterious.values()), test_size)

        fam_data = self.generate_fam_file(control_size, test_size, male_odds, deleterious_group_list)
        main_file = self.population_dir + "population.vcf.gz"

        with bgzf.BgzfWriter(filename=main_file, mode='wt+', compresslevel=compression_level) as f:
            header = gen_vcf_header(fam_data)
            f.write(header)
            print("Outputing VCF lines", flush=True)
            self.write_vcf_snps(fam_data, self.ordered_snps, f)

        print("Finished VCF file output.", flush=True)

    @Timer(text="Finished write_vcf_snps chunk Elapsed time: {:0.4f} seconds", logger=print)
    def write_vcf_snps(self, fam_data, snps, file, header=False):
        processes = []
        # Queue is large since there is possible deadlock scenario if a worker thread gets behind by more than
        # the queue size.
        result_q = Queue(100000)

        # Create a process for each split group
        n_processes = self.num_processes
        if len(snps) < n_processes:
            # Small chunk of work so use 1 process
            n_processes = 1
        start_index = 1
        work_chunks = stripe_list(list(enumerate(snps, start=start_index)), n_processes)

        for i in range(n_processes):
            p = Process(target=self.queue_vcf_snps, args=(fam_data, work_chunks[i], result_q))
            processes.append(p)
            p.start()
        cur_snp = start_index
        backlog = []
        while any(p.is_alive() for p in processes):
            while not result_q.empty():
                snp_tuple = result_q.get()
                if snp_tuple[0] == cur_snp:
                    file.write(snp_tuple[1])
                    cur_snp += 1
                    while backlog and cur_snp == backlog[0][0]:
                        # Next item is on the min-heap
                        heap_tuple = heapq.heappop(backlog)
                        file.write(heap_tuple[1])
                        cur_snp += 1
                else:
                    heapq.heappush(backlog, snp_tuple)
                if snp_tuple[0] % 5000 == 0:
                    print("Output %i/%i VCF lines in file." % (snp_tuple[0], len(snps)), flush=True)

    def queue_vcf_snps(self, fam_data, work_q, result_q):

        num_samples = len(fam_data)
        for snp_num, snp in work_q:
            try:
                #snp_num, snp = work_q.get_nowait()

                sample_values = []
                # Roll the dice for each sample and each allele.
                randoms = numpy.random.rand(num_samples * 2)
                for i, sample in enumerate(fam_data):
                    is_male = sample.is_male()

                    if not is_male and snp.chromosome == 'Y':
                        # No Y chromosome for women
                        sample_values.append(".")
                    if sample.is_control or snp.id not in sample.deleterious_snps:
                        random_roll = randoms[i * 2]
                        selected_nt = snp.pick_allele_index(random_roll)
                        if is_haploid(snp.chromosome, is_male):
                            sample_values.append(str(selected_nt))
                            continue
                        else:
                            random_roll = randoms[i * 2 + 1]
                            other_nt = snp.pick_allele_index(random_roll)
                        sample_values.append("%i/%i" % (selected_nt, other_nt))
                    else:
                        if is_haploid(snp.chromosome, is_male):
                            sample_values.append("1")
                        else:
                            sample_values.append("1/1")
                        # TODO make it so deleterious mutations can be recessive or dominant
                # Output row - CHROM, POS, ID, REF, ALT, QUAL FILTER, INFO, FORMAT, (SAMPLE ID ...)
                # 1      10583 rs58108140  G   A   25   PASS    .    GT     0/0     0/0     0/0
                line = "%s\t%i\trs%s\t%s\t%s\t40\tPASS\t.\tGT\t" % (snp.chromosome,
                                                                    snp.position,
                                                                    snp.id,
                                                                    snp.ref_allele_tuple()[0],
                                                                    snp.alt_alleles()) + \
                       "\t".join(sample_values) + "\n"
                result_q.put((snp_num, line))
            except queue.Empty:
                # pause to allow items being queued to complete being sent
                print("Work Queue empty.")
                time.sleep(1)
                break

    def load_deleterious(self):
        with open(self.deleterious_list_path, 'rt') as f:
            for line in f:
                pg = DeleteriousGroup.from_json(line)
                self.deleterious[pg.name] = pg

    def pick_deleterious_snps(self, snp_data, deleterious_config):
        """
        Pick and store the snps that are the deleterious based on deleterious_config.
        :param deleterious_config: file path for deleterious yaml file
        :param snp_data: SNPTuples which are candidates for being a deleterious
        :return: nothing... self.deleterious is populated
        """

        with open(deleterious_config, 'r') as p:
            deleterious_yml = load(p, Loader=Loader)
            for group, group_attr in deleterious_yml.items():
                iterations = 1
                if group_attr['num_instances']:
                    iterations = int(group_attr['num_instances'])
                for i in range(0, iterations):
                    path_group = DeleteriousGroup.from_yml(group_attr, snp_data, "%s-%s" % (group, i))
                    self.deleterious[path_group.name] = path_group
        with open(self.population_dir + "deleterious.json", 'w') as f:
            for deleterious_group in self.deleterious.values():
                f.write(deleterious_group.to_json() + "\n")


class DeleteriousGroup:

    def __init__(self, name, population_weight):

        self.deleterious = {}
        self.name = name
        self.population_weight = population_weight

    @classmethod
    def init_with_snps(cls, name, mutation_weights, snp_data, population_weight,
                       min_minor_allele_freq=0, max_minor_allele_freq=1.1):
        """
        Inits the deleterious dictionary to be random alleles matching the freq filters. deleterious dict
        stores snp.id => mutation weight mapping.
        :param mutation_weights: list of floats of the value each picked
        :param snp_data: snp dictionary
        :param population_weight: the weight this deleterious group has (shares in test population)
        :param name: name of this group
        :param min_minor_allele_freq: Filter for picking snps. Min MAF of SNP to be in deleterious group
        :param max_minor_allele_freq: Filter for picking snps. Max MAF of SNP to be in deleterious group
        """
        deleterious_grp = cls(name, population_weight)

        filter_snps = min_minor_allele_freq > 0 or max_minor_allele_freq < 0.5
        filtered_list = snp_data
        if filter_snps:
            filtered_list = filter(
                lambda x: min_minor_allele_freq <= (x.minor_allele_tuple()[1]
                                                    - x.ref_allele_tuple()[1]) <= max_minor_allele_freq,
                filtered_list)
        snp_id_list = list(map(lambda x: x.id, filtered_list))
        if len(snp_id_list) == 0:
            raise Exception("All SNPs filtered out. No snps match deleterious filter %f <= freq <= %f" %
                            (min_minor_allele_freq, max_minor_allele_freq))
        i = 0
        for snp_id in numpy.random.choice(a=snp_id_list, size=len(mutation_weights), replace=False):
            deleterious_grp.deleterious[int(snp_id)] = mutation_weights[i]
            i += 1
        return deleterious_grp

    @classmethod
    def from_yml(cls, yml_attr, snp_data, name):
        min_minor_allele_freq = 0
        max_minor_allele_freq = 1
        if yml_attr.get('min_minor_allele_freq'):
            if 0 < yml_attr['min_minor_allele_freq'] < 0.5:
                min_minor_allele_freq = yml_attr['min_minor_allele_freq']
            else:
                raise Exception('min_minor_allele_freq must be between 0 and 0.5. yml value = {}'.format(
                    yml_attr['min_minor_allele_freq']))
        if yml_attr.get('max_minor_allele_freq'):
            if 0 < yml_attr['max_minor_allele_freq'] < 0.5:
                max_minor_allele_freq = yml_attr['max_minor_allele_freq']
            else:
                raise Exception('max_minor_allele_freq must be between 0 and 0.5. yml value = {}'.format(
                    yml_attr['max_minor_allele_freq']))

        return cls.init_with_snps(name,
                                  yml_attr['mutation_weights'],
                                  snp_data,
                                  yml_attr['population_weight'],
                                  min_minor_allele_freq,
                                  max_minor_allele_freq)

    def to_json(self):
        return json.dumps(vars(self))

    @classmethod
    def from_json(cls, json_line):
        pg = json.loads(json_line)
        deleterious_group = cls(pg['name'], pg['population_weight'])
        for snp, weight in pg['deleterious'].items():
            deleterious_group.deleterious[snp] = weight
        return deleterious_group

    def select_mutations(self):
        """
        Randomly select mutations a single individual might have if they are in this DeleteriousGroup
        :return: a colleciton of snp_ids that are mutated
        """
        selected_deleterious = {}
        shuffled_deleterious = list(self.deleterious.items())
        random.shuffle(shuffled_deleterious)  # Shuffle to randomly select
        agg_weight = 0
        for p in shuffled_deleterious:
            selected_deleterious[p[0]] = p[1]
            agg_weight += p[1]  # sum the weights
            if agg_weight >= 1:
                break
        return selected_deleterious


def parse_cmd_args(args):
    arg_parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                         prog='Population Factory',
                                         description='Generates genetic populations using simulated SNP data.')

    arg_parser.add_argument('-s', type=int, dest='size', help='size of afflicted/case group', required=True)
    arg_parser.add_argument('-c', type=int, dest='control_size', help='size of control group', required=True)
    arg_parser.add_argument('-x', type=int, dest='max_snps', help='max number of snps to load/generate')
    arg_parser.add_argument('-p', type=str, default='deleterious.yml',
                            help='location of deleterious config yaml file (default is deleterious.yml)',
                            dest='deleterious_config')
    arg_parser.add_argument('-f', type=float, default=0.005, dest='min_freq',
                            help='min minor allele frequency for a SNP to be included, default is 0.005')

    arg_parser.add_argument("-m", type=float, default=0.5, dest='male_odds',
                            help='odds of a population member being male (default 0.5)')
    arg_parser.add_argument('-n', type=int, default=2, dest='num_processes',
                            help='Number of worker processes to use')
    arg_parser.add_argument('-z', type=int, dest='compression_level', default=6,
                            help='gzip compression level (1=least 9=most) default 6',
                            choices=range(1, 10))
    arg_parser.add_argument('-l', action='store_const', const=False, default=True, dest='generate_snps',
                            help='load from refSNP datababse instead of using '
                                 + 'simulated snps (connection config in db.yml)')
    arg_parser.add_argument('--deleterious_file', type=str,
                            help='<path> to a deleterious.json file that specifies the exact snps to use as deleterious')
    arg_parser.add_argument('--snps_file', type=str,
                            help='<path> location of snps.json.gz file to use as selected snps')
    arg_parser.add_argument('--outdir', type=str,
                            help='<path> directory to use for output files')
    arg_parser.add_argument('--offset', type=int,
                            help='offset to add to all sample ids. Useful for creating VCF files that can be merged')
    return arg_parser.parse_args(args)


def main(sys_args):
    args = parse_cmd_args(sys_args)
    if args.num_processes > 4 and args.compression_level > 5:
        print("Recommend using a compression level of 3 (-z 3) or lower with 5 or more worker threads.")

    if not args.generate_snps:
        db.default_init()
    pop_factory = PopulationFactory(num_processes=args.num_processes,
                                    generate_snps=args.generate_snps,
                                    deleterious_list_path=args.deleterious_file,
                                    sample_id_offset=args.offset,
                                    male_odds=args.male_odds,
                                    deleterious_config=args.deleterious_config,
                                    snps_path=args.snps_file,
                                    output_path=args.outdir)
    pop_factory.generate_population(args.control_size, args.size, args.min_freq,
                                    args.max_snps, args.compression_level)


if __name__ == '__main__':
    main(sys.argv[1:])
