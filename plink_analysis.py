import getopt
import sys


def get_pathogens(pop_path_file):
    """ Returns dict of rsXXXX ids from the pop_pathogens.txt file and the number of cases in the population
    """
    pathogens = {}
    with open(pop_path_file, 'rt') as p:
        for line in p:
            cols = line.split("\t")
            for snp in cols[2:]:
                if snp not in pathogens:
                    pathogens[snp] = 0
                pathogens[snp] += 1
    return pathogens


def output_pathogens_rows(assoc_file, pathogens_dict):
    with open(assoc_file, 'rt') as a:
        i = 0
        pathogen_rows = []
        for line in a:
            i += 1
            if i == 1:
                print(line + "\tCases")
                continue
            cols = line.split()
            p_value = cols[8]
            snp = cols[1]
            if snp in pathogens_dict:
                pathogen_rows.append((line + "\t%i" % pathogens_dict[snp], p_value))
        pathogen_rows.sort(key=lambda x: x[1])
        for row in pathogen_rows:
            print(row[0])


def analyze_assoc_results(argv):
    try:
        opts, args = getopt.getopt(argv, "a:p:", ["help"])
    except getopt.GetoptError as err:
        print(err.msg)
        sys.exit(2)
    plink_assoc = "plink.assoc"
    pathogens_file = "pop_pathogens.txt"
    for opt, arg in opts:
        if opt in "-a":
            plink_assoc = arg
        elif opt in "-p":
            pathogens_file = arg
    pathogens = get_pathogens(pathogens_file)
    # Output pathogens rows. Output Low P values that are not pathogens
    output_pathogens_rows(plink_assoc, pathogens)


if __name__ == '__main__':
    analyze_assoc_results(sys.argv[1:])
