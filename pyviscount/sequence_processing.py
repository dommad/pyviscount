"""Functions related to processing protein/peptide sequences"""
import re



def digest_by_trypsin(sequence, mc):
    """Digestion by trypsin"""

    # cleavage sites and exceptions
    cut_aa = ('K', 'R')
    exception_aa = ('P')

    list_aa = list(sequence)
    cut_sites = [-1]
    length_aa = len(list_aa)

    for idx, item in enumerate(list_aa[:-1]):
        next_aa = list_aa[idx + 1]
        if item in cut_aa and next_aa not in exception_aa:
            cut_sites.append(idx)

    cut_sites.append(length_aa - 1)

    peptides = []
    for start, end in zip(cut_sites, cut_sites[1:]):
        peptides.append(sequence[start + 1 : end + 1])
    
    
    mc_peptides = []
    for cur_mc in range(1, mc + 1):
        mc_peptides.extend(["".join(peptides[i : i + cur_mc + 1]) for i in range(len(peptides) - cur_mc)])

    return peptides + mc_peptides


def clean_sequence(seq):
    """remove unnecessary characters from the peptide sequence"""
    s1 = re.sub(r'{}.*?{}'.format(re.escape("["),re.escape("]")),'', seq)
    s1 = re.sub("-", "", s1)
    if s1[0] == 'n':
        s1 = s1[1:]
    return s1[:-2]


def get_amino_acid_mass_dict():
    """Monoisotopic masses of amino acids"""

    amino_acid_dict = {}
    amino_acid_dict['A'] = 71.03711
    amino_acid_dict['R'] = 156.10111
    amino_acid_dict['N'] = 114.04293
    amino_acid_dict['D'] = 115.02964
    amino_acid_dict['C'] = 103.00919
    amino_acid_dict['F'] = 147.06841
    amino_acid_dict['P'] = 97.05276
    amino_acid_dict['S'] = 87.03203
    amino_acid_dict['T'] = 101.04768
    amino_acid_dict['W'] = 186.07931
    amino_acid_dict['Y'] = 163.06333
    amino_acid_dict['V'] = 99.06841
    amino_acid_dict['E'] = 129.04259
    amino_acid_dict['Q'] = 128.05858
    amino_acid_dict['G'] = 57.02146
    amino_acid_dict['H'] = 137.05891
    amino_acid_dict['I'] = 113.08406
    amino_acid_dict['L'] = 113.08406
    amino_acid_dict['K'] = 128.09496
    amino_acid_dict['M'] = 131.04049
    # amino_acid_dict['X'] = 101.0
    # amino_acid_dict['B'] = 114.5  # C+57
    # amino_acid_dict['J'] = 1  # M+O
    # amino_acid_dict['Z'] = 0.0
    # amino_acid_dict['U'] = 150.95363  # selenocysteine

    return amino_acid_dict