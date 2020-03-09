"""
Common classes used by different python functions
"""

import json
import re

CHROMOSOME_LIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15',
                   '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
CHROMOSOME_PROB = [0.07426087261566,
                   0.07930487311426,
                   0.06669253502772,
                   0.068216704579376,
                   0.060859452377757,
                   0.061620602417568,
                   0.056436996345677,
                   0.052745283940636,
                   0.041811456817423,
                   0.047572674763057,
                   0.046903788666524,
                   0.045558978461098,
                   0.033875108161329,
                   0.030837930905743,
                   0.028329099437382,
                   0.030535626281104,
                   0.026508783521902,
                   0.026711126377244,
                   0.022471493713103,
                   0.021115686613365,
                   0.013429462318399,
                   0.013635819040166,
                   0.048111412615406,
                   0.002454231888101
                   ]

CHROMOSOME_MAX_POSITION = {
    "1": 248946339,
    "2": 242765766,
    "3": 198235509,
    "4": 190181952,
    "5": 181477687,
    "6": 170744571,
    "7": 159335932,
    "8": 145571444,
    "9": 138258771,
    "10": 133787363,
    "11": 135076614,
    "12": 133265032,
    "13": 114352979,
    "14": 107270972,
    "15": 101981181,
    "16": 90228323,
    "17": 83247315,
    "18": 80262386,
    "19": 58607512,
    "20": 64333614,
    "21": 46699955,
    "22": 50806829,
    "X": 156040000,
    "Y": 57217333}


def chromosome_from_filename(filename):
    chr_search = re.search('chr([0-9XYMT]+)', filename, re.IGNORECASE)
    if chr_search:
        return chr_search.group(1)
    else:
        return 'unknown'


def split_list(l, n):
    chunk_size = round(len(l) / n)
    for i in range(n):
        start = i * chunk_size
        if i + 1 == n:
            end = len(l)
        else:
            end = (i + 1) * chunk_size
        yield l[start:end]


def stripe_list(l, num_stripes):
    stripes = []
    for i in range(num_stripes):
        stripes.append([])
    for n, i in enumerate(l):
        stripes[n % num_stripes].append(i)
    return stripes

def obj_from_rowproxy(cls, row_proxy):
    o = cls.__new__(cls)
    for k, v in row_proxy.items():
        try:
            setattr(o, k, v)
        except AttributeError as e:
            print(str(e))
            pass  # Skip unknown attributes
    return o


def is_haploid(chromo, is_male):
    """
    Is this chromosome haploid (one allele per person)
    :param chromo: chromosome letter or number (X, Y, MT or 1-22)
    :param is_male: boolean if male
    :return:
    """
    return (chromo == 'X' and is_male) or chromo == 'MT' or chromo == 'Y'


class Allele:
    """
    Stores data representing an allele of a SNP
    """

    def __init__(self, deleted, inserted, position):
        self.name = Allele.name_string(deleted, inserted)
        self.deleted = deleted
        self.inserted = inserted
        self.position = position
        self.allele_count = 0
        self.total_count = 0
        self.ref_snp_id = None

    def add_observation(self, allele_count, total_count):
        self.allele_count += int(allele_count)
        self.total_count += int(total_count)

    def freq(self):
        return self.allele_count / self.total_count

    def to_dict(self):
        return {"deleted": self.deleted,
                "inserted": self.inserted,
                "position": self.position,
                "allele_count": self.allele_count,
                "total_count": self.total_count}

    @classmethod
    def from_dict(cls, attr_dict):
        a = cls(attr_dict["deleted"], attr_dict["inserted"], attr_dict["position"])
        a.allele_count = attr_dict["allele_count"]
        a.total_count = attr_dict["total_count"]
        return a

    @classmethod
    def from_row_proxy(cls, attr_dict):
        a = cls(attr_dict["deleted"], attr_dict["inserted"], attr_dict["position"])
        a.allele_count = attr_dict["allele_count"]
        return a

    @staticmethod
    def name_string(deleted, inserted):
        return deleted + "->" + inserted


class RefSNP:

    def __init__(self, ref_id, chromosome):
        self.id = ref_id
        self.alleles = []
        self.total_count = None
        self.maf = None
        self.chromosome = chromosome

    def put_allele(self, allele):
        allele.ref_snp_id = self.id
        self.alleles.append(allele)

    def valid_for_plink(self):
        for allele in self.alleles:
            if not allele.inserted or not allele.deleted:
                return False
            if len(allele.deleted) > 1 or len(allele.inserted) > 1:
                # Skip multi-NT snps
                return False
        return True

    def set_maf_and_total_count(self):
        if self.maf:
            return
        self.alleles.sort(key=lambda x: x.allele_count, reverse=True)
        total_count = 0
        for a in self.alleles:
            total_count += a.allele_count
        self.total_count = total_count
        if total_count > 0 and len(self.alleles) > 1:
            self.maf = self.alleles[1].allele_count / total_count

    @classmethod
    def from_json(cls, json_line, chromosome):
        ref_obj = json.loads(json_line)
        ref_snp = cls(ref_obj['id'], str(chromosome))
        for a in ref_obj['alleles']:
            allele = Allele.from_dict(a)
            ref_snp.put_allele(allele)
        ref_snp.set_maf_and_total_count()
        return ref_snp

    @classmethod
    def from_row_proxy(cls, row):
        o = cls(row['id'], row['chromosome'])
        o.total_count = row['total_count']
        o.maf = row['maf']
        return o

    @classmethod
    def from_nih_json(cls, json_line, chromosome):
        ref_obj = json.loads(json_line)
        ref_snp = cls(ref_obj['refsnp_id'], chromosome)
        if 'primary_snapshot_data' in ref_obj:
            placements = ref_obj['primary_snapshot_data']['placements_with_allele']

            for alleleinfo in placements:
                placement_annot = alleleinfo['placement_annot']
                if alleleinfo['is_ptlp'] and \
                        len(placement_annot['seq_id_traits_by_assembly']) > 0:
                    ref_snp.assembly_name = placement_annot[
                        'seq_id_traits_by_assembly'][0]['assembly_name']

                    for a in alleleinfo['alleles']:
                        spdi = a['allele']['spdi']
                        allele = Allele(spdi['deleted_sequence'],
                                        spdi['inserted_sequence'],
                                        spdi['position'])
                        ref_snp.put_allele(allele)
            for allele_annotation in ref_obj['primary_snapshot_data']['allele_annotations']:
                for freq in allele_annotation['frequency']:
                    obs = freq['observation']
                    name = Allele.name_string(obs['deleted_sequence'],
                                              obs['inserted_sequence'])
                    for allele in ref_snp.alleles:
                        if name == allele.name:
                            allele.add_observation(freq['allele_count'], freq['total_count'])
        ref_snp.set_maf_and_total_count()
        return ref_snp

    @classmethod
    def update_total_counts(cls, session):
        """
        Slow function that updates all the total_counts based on sum of all allele counts for this ref_snp
        :param session: db session
        :return:
        """
        update_total_sql = """
        update ref_snps set total_count = 
        (select total_count from 
            (select ref_snp_id, sum(allele_count) as total_count from alleles 
             join ref_snps on ref_snps.id = alleles.ref_snp_id group by ref_snp_id)
         where
         ref_snps.id = ref_snp_id);
        """
        session.execute(update_total_sql)
        session.commit()

    @classmethod
    def update_maf(cls, session):
        update_maf_sql = """
        --- MAF query
        update ref_snps set maf = 
        ( select sec_high * 1.0 / ref_snps.total_count  from 
            ( select a.ref_snp_id, max(a.allele_count) as sec_high from 
                (select alleles.id, ref_snp_id, max(allele_count) as highest 
                 from alleles 
                 join ref_snps on ref_snps.id = alleles.ref_snp_id
                 group by ref_snp_id) as x
              join alleles a on x.ref_snp_id = a.ref_snp_id and a.id != x.id 
              group by a.ref_snp_id) as s
          where ref_snps.id = s.ref_snp_id
        ) ;
        """
        session.execute(update_maf_sql)
        session.commit()

    @classmethod
    def delete_chromosomes(cls, chromosomes, session):
        """
        Delete from DB all allele and refSNP records for the passed in chromosome list
        :param chromosomes: list of chromosomes (strings) to delete
        :param session: db session/connection
        :return: nothing
        """
        chromo_in_clause = ",".join(map(lambda x: "\'%s\'" % x, chromosomes))
        delete_alleles_sql = """
        delete from alleles where ref_snp_id in (select r.id from ref_snps r where chromosome in (%s))
        """ % chromo_in_clause
        session.execute(delete_alleles_sql)

        delete_ref_snp_sql = """
        delete from ref_snps where chromosome in (%s)
        """ % chromo_in_clause
        session.execute(delete_ref_snp_sql)

    def __str__(self):
        json_hash = {"id": self.id}
        if len(self.alleles) > 0:
            json_hash["alleles"] = []
        for allele in self.alleles:
            json_hash["alleles"].append(allele.to_dict)
        return json.dumps(json_hash)
