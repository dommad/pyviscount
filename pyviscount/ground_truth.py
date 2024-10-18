"""Module for processing of validation by partition workflow"""
from xml.etree import ElementTree as ET
from functools import partial
from collections import deque, namedtuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyteomics import pepxml, mzxml
from pyopenms import MSExperiment, MzMLFile, MzXMLFile
import parsers as sp
from Bio import SeqIO


def execute_main():
    """Execute the main script - generate disjoint peptide tables
    for SpectraST to enable search compatible with validation by partition"""

    sptxt = sp.SpTXT()
    sptxt.read_sptxt('/home/dominik/spectral_libraries/chinese_hamster_hcd.sptxt')
    peptides = sptxt.read_only_peptide()

    np.random.shuffle(peptides)
    peptides_1 = peptides[:int(len(peptides)/2)]
    peptides_2 = peptides[int(len(peptides)/2):]

    export_peps(peptides_1, 'hamster_peps1')
    export_peps(peptides_2, 'hamster_peps2')


def digest_fasta_database(inputfile):
    """Digest a .fasta database and """
    handle = SeqIO.parse(inputfile, 'fasta')
    core = inputfile.split('.')[0]
    db_1 = open(core + "_setAB" + '.fasta', 'w', encoding='utf-8')
    # db_2 = open(core + "_setB" + '.fasta', 'w', encoding='uft-8')


    peptides = []
    for record in handle:
        proseq = str(record.seq)
        peptide_list = digest_by_trypsin(proseq, 0)

        for pep in peptide_list:
            peptides.append(pep)

    pep_set = set(peptides)

    idx = 0
    for pep in pep_set:
        toss = np.random.randint(0, 2)
        if toss == 1:
            db_1.write(">" + "setA_" + str(idx) + '\n' + pep + '\n')
        else:
            db_1.write(">" + "setB_" + str(idx) + '\n' + pep + '\n')

        idx += 1

    db_1.close()


def digest_by_trypsin(proseq, miss_cleavage):
    """digest the protein according to trypsin rules"""
    peptides = []
    cut_sites = [0]
    for i in range(0, len(proseq) - 1):
        if proseq[i] == 'K' and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)
        elif proseq[i] == 'R' and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)

    if cut_sites[-1] != len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites) > 2:
        if miss_cleavage == 0:
            for j in range(0, len(cut_sites) - 1):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])

        elif miss_cleavage == 1:
            for j in range(0, len(cut_sites) - 2):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])

            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])

        elif miss_cleavage == 2:
            for j in range(0, len(cut_sites) - 3):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 3]])

            peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
            peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    else:
        peptides.append(proseq)
    return peptides

# def digest_by_trypsin(proseq, miss_cleavage):
#     """Digest the protein according to trypsin rules"""
#     cut_sites = [i + 1 for i in range(len(proseq) - 1) if (proseq[i] in 'KR' and proseq[i + 1] != 'P')] + [len(proseq)]
#     peptides = []
    
#     if miss_cleavage == 0:
#         peptides = [proseq[cut_sites[j]:cut_sites[j + 1]] for j in range(len(cut_sites) - 1)]
#     elif miss_cleavage == 1:
#         for j in range(len(cut_sites) - 1):
#             peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
#             if j < len(cut_sites) - 2:
#                 peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])
#     elif miss_cleavage == 2:
#         for j in range(len(cut_sites) - 1):
#             peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
#             if j < len(cut_sites) - 2:
#                 peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])
#             if j < len(cut_sites) - 3:
#                 peptides.append(proseq[cut_sites[j]:cut_sites[j + 3])
    
#     return peptides

def digest_database(inputfile):
    """digest the fasta database"""
    handle = SeqIO.parse(inputfile, 'fasta')
    core = inputfile.split('.')[0]
    output = open(core + '.csv', 'w')

    for record in handle:
        proseq = str(record.seq)
        peptide_list = digest_by_trypsin(proseq, 2)
        for peptide in peptide_list:
            if len(peptide) >= 7:
                output.write(record.id + '\t' + peptide + '\n')



def get_pep_mass_dict(min_iprob, filename, engine='SpectraST'):

    """Generate a map of peptides to their masses

    Parameters
    ----------
    min_iprob : float
        minimum probability threshold
    filename : str
        name of the interact-pep.xml file to be processed

    Returns
    -------
    spectra: dict
        a dictionary of peptide/charge (SpectraST format): precursor neutral mass
    """

    file = pepxml.read(filename)
    spectra = {}

    for spec in file:

        prob = spec['search_hit'][0]['analysis_result'][0]['peptideprophet_result']['probability']
        if prob < min_iprob:
            continue

        #scan = spec["start_scan"]
        charge = spec["assumed_charge"]
        mass = spec["precursor_neutral_mass"]
        peptide = spec["search_hit"][0]["modified_peptide"]
        key = f"{peptide}/{charge}"
        if engine == 'Comet':
            key = spec["search_hit"][0]["peptide"]

        if spectra.get(key) is None:
            spectra[key] = mass

    return spectra

def plot_pep_mass_distr(spectra, pep_dict):
    """plot mass distribution of spectra to verify the split is not biased"""

    peps = np.array(tuple(spectra.keys()))
    np.random.shuffle(peps)
    pep_a = peps[:int(len(peps)/2)]
    pep_b = peps[int(len(peps)/2):]

    mass1 = list(map(lambda x: pep_dict[x], pep_a))
    mass2 = list(map(lambda x: pep_dict[x], pep_b))

    plt.hist(mass1, density=True, bins=100, alpha=0.5)
    plt.hist(mass2, density=True, bins=100, alpha=0.5)
    #plt.hist(masses[sample],  density=True, bins=100)
    plt.show()
    plt.savefig("./mass_distr.png", dpi=500, bbox_inches="tight")



def split_peptides(spectra):

    """Randomly split the set of peptides into set A and set B

    Parameters
    ----------
    spectra : dict
        a dictionary of peptide/charge (SpectraST format): precursor neutral mass

    Returns
    -------
    (list, list)
        tuple of two peptide lists for set A and set B, respectively
    """

    peps = np.array(tuple(spectra.keys()))
    np.random.shuffle(peps)
    peps1 = peps[:int(len(peps)/2)]
    peps2 = peps[int(len(peps)/2):]

    return (peps1, peps2)

def export_peps(peps, outname):

    """Export peptide list that can be fed to SpectraST

    Parameters
    ----------
    peps : list
        list of peptides in SpectraST-compatible format

    Returns
    -------
    
    """

    dfr = pd.DataFrame(peps)
    dfr.to_csv(f"{outname}.tsv", sep='\t', header=None, index=None)



def make_ground_truth(filename, min_iprob, pep_dict):

    """Get dictionary of spectra scan: peptide (SpectraST format)
    only for IDs that belong to either set A or set B from spectral library

    Parameters
    ----------
    min_iprob : float
        minimum probability threshold for PSMs
    filename : str
        name of the interact-pep.xml file to be processed
    pep_dict : dict
        dictionary of peptide set A and B from spectral library

    Returns
    -------
    dict
        dictionary of scan: peptide (SpectraST format)
    """

    scan_pep_dict = {}
    file = pepxml.read(filename)
    ref_peps = pep_dict.keys()

    for spec in file:

        prob = spec['search_hit'][0]['analysis_result'][0]['peptideprophet_result']['probability']
        if prob < min_iprob:
            continue

        scan = spec["start_scan"]
        charge = spec["assumed_charge"]
        peptide = spec["search_hit"][0]["modified_peptide"]
        value = f"{peptide}/{charge}"

        if value not in ref_peps:
            continue
        if scan_pep_dict.get(scan) is None:
            scan_pep_dict[scan] = value

    return scan_pep_dict


def get_filtered_mgf(gt_dict, input_file, output_name):

    """"input: mzXML file, gt dictionary

    need to connect scan number with ID in mzXML file
    """
    gt_keys = gt_dict.keys()
    mzxml_file = mzxml.MzXML(input_file)
    for spectrum in mzxml_file:

        with open(f"{output_name}.mgf", 'a') as the_file:
            if spectrum["msLevel"] == 2:
                scan = int(spectrum["id"])

                if scan not in gt_keys:
                    continue

                ret_time = float(spectrum["retentionTime"].real)
                pep_mz = float(spectrum["precursorMz"][0]["precursorMz"])
                charge = spectrum["precursorMz"][0]["precursorCharge"]
                mzs = spectrum["m/z array"]
                ints = spectrum["intensity array"]

                the_file.write("BEGIN IONS\n")
                the_file.write(f"TITLE={output_name}.{scan}.{scan}.{charge}\n")
                the_file.write(f"SCAN={scan}\n")
                the_file.write(f"RTINSECONDS={ret_time:.4f}\n")
                the_file.write(f"PEPMASS={pep_mz:.4f}\n")
                the_file.write(f"CHARGE={charge}+\n")

                for idx, mz_val in mzs:
                    the_file.write(f"{mz_val:.4f} {ints[idx]:.4f}\n")

                the_file.write("END IONS\n")


def compare_gt(gt_dict, results_file, min_prob):

    """Compare SpectraST search results against the ground truth dict

    Parameters
    ----------
    gt_dict : dict
        ground truth dictionary (scan: peptide)
    results_file : str
        name of the SpectraST pepxml file with the results to be processed
    min_prob : float
        minimum probability of iProphet/PeptideProphet to filter the PSMs

    Returns
    -------
    dict
        dictionary of scan: identification label
    (0: incorrect ID, 1: correct ID, 2: ID outside the ground truth)
    """

    spectrast = pepxml.read(results_file)
    results = {}

    for spec in spectrast:

        prob = spec['search_hit'][0]['analysis_result'][0]['peptideprophet_result']['probability']
        if prob < min_prob:
            continue

        scan = spec["start_scan"]
        charge = spec["assumed_charge"]
        peptide = spec["search_hit"][0]["modified_peptide"]
        value = f"{peptide}/{charge}"
        # print(f"this is value: {value}")

        if gt_dict.get(scan) is None:
            results[scan] = 2
        else:
            if gt_dict[scan] == value:
                results[scan] = 1
            else:
                results[scan] = 0

        # print(f"this is gt: {results[scan]}")

    return results, spectrast


Psm = namedtuple("Psm", ["spec_name",
                         "prob", "label", "score",
                         "exp_pep", "gt_pep", "decoy"])

def peptideprophet_validation(interact_file, gt_dict):
    """Extract prob values and FDR thresholds from PeptideProphet results

    Parameters
    ----------
    interact_file : str
        path to interact.pep.xml file with PeptideProphet/iProphet results
    gt_dict : dict
        ground truth dictionary, scan: peptide

    Returns
    -------
    dfs : DataFrame
        DataFrame constructed from Psm namedtuple records
    thr : list
        PeptideProphet probability values for selected FDR thresholds
    """

    records = deque()
    spec_data = pepxml.read(interact_file)

    for res in spec_data:
        if 'search_hit' in res.keys():
            decoy_status = bool('DECOY' in res['search_hit'][0]['proteins'][0]['protein'])
            top = res['search_hit'][0]
            exp_pep = f"{top['modified_peptide']}/{res['assumed_charge']}"
            gt_pep = gt_dict[int(res['start_scan'])]

            cur_psm = Psm(res['spectrum'],
                          top['analysis_result'][0]['peptideprophet_result']['probability'],
                          bool(exp_pep == gt_pep),
                          top['search_score']['dot'],
                          exp_pep,
                          gt_pep,
                          decoy_status)

            records.append(cur_psm)

    dfs = pd.DataFrame(records)
    #dfs = dfs[dfs.decoy == 0]
    dfs = dfs.sort_values('prob', inplace=False, ascending=True)
    dfs = dfs.reset_index(drop=True)
    dfs.index += 1
    tree = ET.parse(interact_file)
    root = tree.getroot()
    no_files = 1
    # 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05
    fdr_idx = [34, 40, 45, 46, 47, 48, 49, 50, 51]
    thr = list(map(lambda x: float(root[0][0][int(no_files)][x].attrib['min_prob']), fdr_idx))
    return dfs, thr


def get_fdp(threshold, dfr):
    """Calculates FDP based on PeptideProphet probability values"""
    discovered = dfr[dfr['prob'] >= threshold]
    incorr = len(discovered[discovered.label == False])
    return incorr/len(discovered)


def plot_fdp_fdr(dfr, thrs, outname):
    """Plots FDP vs. FDR thresholds"""

    fdps = list(map(partial(get_fdp, dfr=dfr), thrs))
    fig, axes = plt.subplots(figsize=(5,5))
    fdrs = [0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05]
    axes.scatter(fdrs, fdps, s=5)
    axes.plot(fdrs, fdps)
    axes.plot([0, 0.06], [0, 0.06], c='grey')
    axes.set_xlim(0, 0.05)
    axes.set_ylim(0, max(fdps)+0.005)
    axes.set_xlabel("FDR threshold")
    axes.set_ylabel("FDP")
    fig.savefig(f'{outname}.png', dpi=500, bbox_inches='tight')


def get_filtered_mzml(gt_dict, input_file, output_file):
    """Create mzML with spectra present in ground truth dict"""

    gt_keys = gt_dict.keys()
    mzml_dict = MSExperiment()
    if input_file.split('.')[-1]  == 'mzML':
        MzMLFile().load(input_file, mzml_dict)
        ext = "mzML"
    if input_file.split('.')[-1]  == 'mzXML':
        MzXMLFile().load(input_file, mzml_dict)
        ext = "mzXML"

    new_mzml = MSExperiment()
    in_gt = filter(partial(get_conditions_mzml, gt_keys=gt_keys), mzml_dict.getSpectra())

    for spec in in_gt:
        new_mzml.addSpectrum(spec)

    if ext == 'mzXML':
        MzXMLFile().store(f'filtered_{output_file}.{ext}', new_mzml)
    else:
        MzMLFile().store(f'filtered_{output_file}.{ext}', new_mzml)


def get_conditions_mzml(spec, gt_keys):
    """Conditions for filtering out MS2 spectra found in GT dict"""
    is_ms2 = bool(spec.getMSLevel() == 2)
    is_in_gt = bool(int(str(spec.getNativeID()).split('=')[-1]) in gt_keys)
    return is_ms2 and is_in_gt
