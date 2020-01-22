"""
Common classes used by different python functions
"""

import json

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
    refsnp -> allele
    id        (id, refsnp_id, deleted, inserted, seq_id, position, allele_count)
    chromo
    MAF
    total_count
    """

    def __init__(self, deleted, inserted, seq_id, position):
        self.name = Allele.name_string(deleted, inserted)
        self.deleted = deleted
        self.inserted = inserted
        self.seq_id = seq_id
        self.position = position
        self.allele_count = 0
        self.total_count = 0

    def add_observation(self, allele_count, total_count):
        self.allele_count += int(allele_count)
        self.total_count += int(total_count)

    def freq(self):
        return self.allele_count / self.total_count

    def to_dict(self):
        return {"deleted": self.deleted,
                "inserted": self.inserted,
                "position":  self.position,
                "seq_id": self.seq_id,
                "allele_count": self.allele_count,
                "total_count": + self.total_count}

    @staticmethod
    def name_string(deleted, inserted):
        return deleted + "->" + inserted


class RefSNP:

    def __init__(self, ref_id):
        self.id = ref_id
        self.alleles = {}

    def put_allele(self, allele):
        self.alleles[allele.name] = allele

    @classmethod
    def from_json(cls, json_line):
        ref_obj = json.loads(json_line)
        ref_snp = cls(ref_obj['refsnp_id'])
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
                                        spdi['seq_id'],
                                        spdi['position'])
                        ref_snp.put_allele(allele)
            for allele_annotation in ref_obj['primary_snapshot_data']['allele_annotations']:
                for freq in allele_annotation['frequency']:
                    obs = freq['observation']
                    name = Allele.name_string(obs['deleted_sequence'],
                                              obs['inserted_sequence'])
                    if name in ref_snp.alleles:
                        ref_snp.alleles[name].add_observation(freq['allele_count'], freq['total_count'])
        return ref_snp

    def total_allele_count(self):
        sum_count = 0
        for a in self.alleles.values():
            sum_count += a.allele_count
        return sum_count

    def __str__(self):
        json_hash = {"id": self.id}
        if len(self.alleles) > 0:
            json_hash["alleles"] = []
        for allele in self.alleles.values():
            json_hash["alleles"].append(allele.to_dict)
        return json.dumps(json_hash)
