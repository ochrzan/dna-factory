"""
Common classes used by different python functions
"""

import json
import re


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
                "seq_id": self.seq_id,
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

    @classmethod
    def from_json(cls, json_line, chromosome):
        ref_obj = json.loads(json_line)
        ref_snp = cls(ref_obj['id'], str(chromosome))
        for a in ref_obj['alleles']:
            allele = Allele.from_dict(a)
            ref_snp.put_allele(allele)
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
        return ref_snp

    @classmethod
    def update_total_counts(cls, chromosome, session):
        """
        Slow function that updates all the total_counts based on sum of all allele counts for this ref_snp
        :param chromosome: The chromosome data to update
        :param session: db session
        :return:
        """
        update_total_sql = """
        update ref_snps set total_count = 
        (select total_count from 
            (select ref_snp_id, sum(allele_count) as total_count from alleles 
             join ref_snps on ref_snps.id = alleles.ref_snp_id and ref_snps.chromosome = '%s' group by ref_snp_id)
         where
         ref_snps.id = ref_snp_id);
        """ % chromosome
        session.execute(update_total_sql)
        session.commit()

    @classmethod
    def update_maf(cls, chromosome, session):
        update_maf_sql = """
        --- MAF query
        update ref_snps set maf = 
        ( select sec_high * 1.0 / ref_snps.total_count  from 
            ( select a.ref_snp_id, max(a.allele_count) as sec_high from 
                (select alleles.id, ref_snp_id, max(allele_count) as highest 
                 from alleles 
                 join ref_snps on ref_snps.id = alleles.ref_snp_id and ref_snps.chromosome = '%s'
                 group by ref_snp_id) as x
              join alleles a on x.ref_snp_id = a.ref_snp_id and a.id != x.id 
              group by a.ref_snp_id) as s
          where ref_snps.id = s.ref_snp_id
        ) ;
        """ % chromosome
        session.execute(update_maf_sql)
        session.commit()

    def total_allele_count(self):
        sum_count = 0
        for a in self.alleles:
            sum_count += a.allele_count
        return sum_count

    def __str__(self):
        json_hash = {"id": self.id}
        if len(self.alleles) > 0:
            json_hash["alleles"] = []
        for allele in self.alleles:
            json_hash["alleles"].append(allele.to_dict)
        return json.dumps(json_hash)



