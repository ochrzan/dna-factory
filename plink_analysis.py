import getopt
import sys


def get_deleterious(pop_path_file):
    """ Returns dict of rsXXXX ids from the pop_deleterious.txt file and the number of cases in the population
    """
    deleterious = {}
    with open(pop_path_file, 'rt') as p:
        for line in p:
            cols = line.split("\t")
            for snp in cols[2:]:
                snp_key = snp.strip()
                if snp_key not in deleterious:
                    deleterious[snp_key] = 0
                deleterious[snp_key] += 1
    return deleterious


def output_deleterious_rows(assoc_file, deleterious_dict):
    print('******* Plink assoc lines for SNPs that are deleterious *****')
    with open(assoc_file, 'rt') as a:
        i = 0
        deleterious_rows = []
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
            if snp in deleterious_dict:
                deleterious_rows.append((line.replace("\n", "") + "\t%i" % deleterious_dict[snp], p_value))
        deleterious_rows.sort(key=lambda x: x[1])
        for row in deleterious_rows:
            print(row[0])


def output_low_p_vals(assoc_file, deleterious_dict, num_vals=30):
    print('******* SNPs with lowest P Values. SNPs with an "*" are deleterious *****')
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
            if snp in deleterious_dict:
                cols.append(str(deleterious_dict[snp]))
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
    # TODO switch to argparse
    try:
        opts, args = getopt.getopt(argv, "a:p:", ["help"])
    except getopt.GetoptError as err:
        print(err.msg)
        sys.exit(2)
    plink_assoc = "plink.assoc"
    deleterious_file = "pop_deleterious.txt"
    for opt, arg in opts:
        if opt in "-a":
            plink_assoc = arg
        elif opt in "-p":
            deleterious_file = arg
    deleterious = get_deleterious(deleterious_file)
    # Output deleterious rows. Output Low P values that are not deleterious
    output_deleterious_rows(plink_assoc, deleterious)
    output_low_p_vals(plink_assoc, deleterious)


if __name__ == '__main__':
    analyze_assoc_results(sys.argv[1:])
