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
                snp_key = snp.strip()
                if snp_key not in pathogens:
                    pathogens[snp_key] = 0
                pathogens[snp_key] += 1
    return pathogens


def output_pathogens_rows(assoc_file, pathogens_dict):
    print('******* Plink assoc lines for SNPs that are pathogens *****')
    with open(assoc_file, 'rt') as a:
        i = 0
        pathogen_rows = []
        for line in a:
            i += 1
            if i == 1:
                cols = line.split()
                for j, label in enumerate(cols):
                    if label == "P":
                        p_index = j
                    if label == "ID":
                        id_index = j
                print(line.replace("\n", "") + "\tCases")
                continue
            cols = line.split()
            p_value = cols[p_index]
            snp = cols[id_index]
            if snp in pathogens_dict:
                pathogen_rows.append((line.replace("\n", "") + "\t%i" % pathogens_dict[snp], p_value))
        pathogen_rows.sort(key=lambda x: x[1])
        for row in pathogen_rows:
            print(row[0])


def output_low_p_vals(assoc_file, pathogens_dict, num_vals=30):
    print('******* SNPs with lowest P Values. SNPs with an "*" are pathogens *****')
    with open(assoc_file, 'rt') as a:
        i = 0
        all_rows = []
        for line in a:
            i += 1
            if i == 1:
                print(line.replace("\n", "") + "\tCases")
                cols = line.split()
                for j, label in enumerate(cols):
                    if label == "P":
                        p_index = j
                    if label == "ID":
                        id_index = j
                continue
            cols = line.split()
            snp = cols[id_index]
            if snp in pathogens_dict:
                cols.append(str(pathogens_dict[snp]))
                cols.append("*PATHOGEN*")
            else:
                cols.extend(("", ""))
            all_rows.append(cols)
        all_rows.sort(key=lambda x: x[p_index])
        for i, row in enumerate(all_rows):
            if i >= num_vals:
                break
            print("\t".join(row))


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
    output_low_p_vals(plink_assoc, pathogens)


if __name__ == '__main__':
    analyze_assoc_results(sys.argv[1:])
