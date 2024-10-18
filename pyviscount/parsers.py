"""Module for processing data from spectral libraries in format .sptxt"""
import re
from typing import List
from abc import ABC, abstractmethod
from dataclasses import dataclass
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from .constants import TH_BETA, TH_N0
from .utils import _is_numeric, ParserError

FILE_FORMATS = ['pep.xml', 'txt', 'mzid', 'pepxml']
SEARCH_ENGINES = ['Comet', 'SpectraST', 'Tide', 'MSFRagger', 'MSGF+']


class PSMParser(ABC):
    def __init__(self, decoy_tag='decoy'):
        self.decoy_tag = decoy_tag


    @abstractmethod
    def rename_columns(self):
        pass


    def parse(self, file_name):
        #self.decoy_tag = config.get("validation.extra", "decoy_tag", fallback='decoy').strip()
        after_dots = file_name.lower().split('/')[-1].split('.')
        
        if after_dots[-2:] == ['pep', 'xml']:
            file_ext = 'pepxml'
        elif after_dots[-1] in FILE_FORMATS:
            file_ext = after_dots[-1].lower()
        else:
            raise ValueError(f"Unsupported file extension for the file: {file_name}")

        return getattr(self, f"parse_{file_ext}")(file_name)


    def parse_pepxml(self, file_name):
        """Parses pepxml (Comet) and outputs pandas dataframe"""
        #TODO: hit numbers not updating for Tide pep.xml
        try:

            # Define the namespace used in the XML
            ns = {'pepXML': 'http://regis-web.systemsbiology.net/pepXML'}

            def get_spectrum_queries():
                for event, elem in ET.iterparse(file_name, events=('start', 'end')):
                    if event == 'start' and elem.tag.endswith('spectrum_query'):
                        spectrum_info = elem.attrib
                        spectrum_info = self.turn_strings_into_floats(spectrum_info)

                    elif event == 'end' and elem.tag.endswith('search_hit'):
                        search_hit_info = self.parse_spectrum_info(elem.iter())
                        analysis_result_list = elem.findall('.//pepXML:analysis_result', namespaces=ns)
                        analysis_result_info = {}

                        for analysis_result in analysis_result_list:
                            analysis_tag = analysis_result.attrib.get('analysis', '')
                            raw_analysis_info = self.parse_spectrum_info(analysis_result.iter())
                            cur_analysis = {f"{analysis_tag}_{key}": val for key, val in raw_analysis_info.items()}
                            analysis_result_info.update(cur_analysis)

                        combined = {**spectrum_info, **search_hit_info, **analysis_result_info}
                        yield combined
                        elem.clear()

            df = pd.DataFrame(get_spectrum_queries())
            df.rename(columns=self.rename_columns(), inplace=True)
            df = df.fillna(0)
            df['modifications'] = list(zip(df['position'], df['mass']))

            return self.add_extra_columns(df)
        
        except ParserError as err:
            print(f"Error parsing file: {err}")
            return None


    def turn_strings_into_floats(self, x):
        return {k: float(v) if _is_numeric(v) else v for k, v in x.items()}
        

    def parse_spectrum_info(self, spectrum_info):
        master_dict = {}
        for item in spectrum_info:
            cur_dict = item.attrib
            #if all(key in cur_dict for key in ['name', 'value']):
            if 'name' in cur_dict and 'value' in cur_dict:
                master_dict[cur_dict['name']] = cur_dict['value']
                del cur_dict['name']
                del cur_dict['value']

            master_dict.update(cur_dict)

        return self.turn_strings_into_floats(master_dict)


    def parse_txt(self, file_name, sep='\t'):
        # Implement common parsing logic
        df = pd.read_csv(file_name, sep=sep)
        # renaming will be handled by individual engines
        df.rename(columns=self.rename_columns(), inplace=True)
        df = df.fillna(0)

        return self.add_extra_columns(df)


    def parse_mzid(self, xml_file_path):
      
        def get_spectra():
            peptides = {}
            peptide_evidence = {}
            db_proteins = {}
            # Iterate over XML elements as they are parsed
            for event, elem in ET.iterparse(xml_file_path, events=('start', 'end')):
                if event == 'start':
                    
                    if elem.tag.endswith("MzIdentML"):
                        ns = {'mzIdentML': elem.tag.split('}')[0][1:]}
                    
                    elif elem.tag.endswith("DBSequence"):
                        protein_info = elem.attrib
                        protein_info['protein_name'] = protein_info['accession']
                        prot_dict = {protein_info['id']: protein_info}
                        db_proteins.update(prot_dict)

                    # Process start events to collect Peptide information
                    elif elem.tag.endswith("Peptide") and elem.attrib.get("id"):
                        pep_info = self.parse_spectrum_info(elem.iter())
                        pep_dict = {pep_info['id']: pep_info}
                        peptides.update(pep_dict)
                    
                    elif elem.tag.endswith("PeptideEvidence") and elem.attrib.get("id"):
                        evidence_dict = self.parse_spectrum_info(elem.iter())
                        peptide_evidence.update({evidence_dict['id']: evidence_dict})
                        

                elif event == 'end':
                    # Process end events to collect SpectrumIdentificationResults
                    combined_info = {}
                    if elem.tag.endswith("SpectrumIdentificationResult") and elem.attrib.get("id"):
                        spectrum_attrib = elem.attrib

                        for spectrum_identification_item in elem.findall('.//mzIdentML:SpectrumIdentificationItem', ns):
                            psm_info = self.parse_spectrum_info(spectrum_identification_item.iter())
                            mod_info = peptides.get(psm_info['peptide_ref'], {})
                            # TODO: consider moving the logic of protein info to an abstract function, Comet handles it differently than MSGF+
                            peptide_evidence_info = peptide_evidence.get(psm_info.get('peptideEvidence_ref', ""), "")
                            protein_info = db_proteins.get(peptide_evidence_info.get('dBSequence_ref', ""), {})
                            combined_info = {**spectrum_attrib, **psm_info, **mod_info, **protein_info, **peptide_evidence_info}
                            yield combined_info

                        elem.clear()

        # Create a DataFrame
        df = pd.DataFrame(get_spectra())
        df.rename(columns=self.rename_columns(), inplace=True)

        return self.add_extra_columns(df)


    def add_extra_columns(self, df):

        df.loc[:, 'is_decoy'] = df['protein'].str.lower().str.contains(self.decoy_tag)
        df.loc[:, 'tev'] = self.calculate_tev(df, -TH_BETA, TH_N0)

        return df


    @staticmethod
    def calculate_tev(df: pd.DataFrame, par_a: float, par_n0: float) -> pd.Series:
        """
        Calculate the log-transformed e-value (TEV) score based on the given parameters.

        Parameters:
        - df (pd.DataFrame): Input DataFrame containing relevant information.
        - par_a (float): The 'a' parameter used in TEV score calculation.
        - par_n0 (float): The 'N0' parameter used in TEV score calculation.

        Returns:
        np.ndarray: An array containing TEV scores for each row in the DataFrame.
        """

        if 'e_value' in df.columns:
            e_values = np.maximum(df['e_value'].values, 1e-16)
            return par_a * np.log(e_values / par_n0)
        
        df['p_value'].replace(-1, 1000, inplace=True)
        p_values = np.maximum(df['p_value'].values, 1e-16)
        num_candidates = np.maximum(df['num_candidates'].values, 1)
        return par_a * np.log(p_values * num_candidates / par_n0)



class CometParser(PSMParser):
    def __init__(self):
        super().__init__()

    def rename_columns(self):
        # Define column renaming logic specific to Comet
        columns = { # pepxml columns
                    'start_scan': 'scan',
                    'peptide': 'sequence',
                    'num_matched_peptides': 'num_candidates',
                    'expect': 'e_value',
                    'assumed_charge': 'charge',
                    'modified_peptide': 'modifications',
                    #txt columns
                    'e-value': 'e_value',
                    'plain_peptide': 'sequence',
                    # mzid columns
                    'spectrumID': 'scan',
                    'rank': 'hit_rank',
                    'chargeState': 'charge',
                    'peptide_ref': 'sequence',
                    'Comet:expectation value': 'e_value',
                    'peptideEvidence_ref': 'protein',
                    'Comet:xcorr': 'xcorr',
                    }
        return columns


class SpectraSTParser(PSMParser):
    def __init__(self):
        super().__init__()

    def rename_columns(self):
        # Define column renaming logic specific to SpectraST
        columns = {'start_scan': 'scan',
                    'peptide': 'sequence',
                    'hits_num': 'num_candidates',
                    }
        return columns
    

class TideParser(PSMParser):
    def __init__(self):
        super().__init__()

    def rename_columns(self):
        # Define column renaming logic specific to Tide
        columns = {'exact p-value': 'p_value',
                    'distinct matches/spectrum': 'num_candidates',
                    'xcorr rank': 'hit_rank',
                    'protein id': 'protein',
                    'refactored xcorr': 'refactored_xcorr',
                    
                    #pepxml columns
                    'exact_pvalue': 'p_value',
                    'num_matched_peptides': 'num_candidates',
                    'peptide': 'sequence'}
        
        return columns


class MSFraggerParser(PSMParser):
    def __init__(self):
        super().__init__()

    def rename_columns(self):
        # Define column renaming logic specific to MSFragger
        columns = {'SpectrumID': 'scan',
                'Rank': 'hit_rank', 
                'Peptide_Sequence': 'sequence',
                'Modifications': 'modifications',
                'expect': 'e_value',
                'start_scan': 'scan',
                'peptide': 'sequence'
                }

        return columns
    

class MSGFParser(PSMParser):
    def __init__(self):
        super().__init__()

    def rename_columns(self):
        # Define column renaming logic specific to MSGF+
        columns = {'start': 'scan',
                    'peptide_ref': 'sequence',
                    'protein_name': 'protein',
                    'hits_num': 'num_candidates',
                    'MS-GF:EValue': 'e_value',
                    'monoisotopicMassDelta': 'modifications',
                    'MS-GF:RawScore': 'msgf_raw_score',
                    'MS-GF:DeNovoScore': 'msgf_denovo_score',
                    'MS-GF:SpecEValue': 'msgf_spec_e_value',
                    'chargeState': 'charge',
                    'rank': 'hit_rank'
                    }
        return columns


class ParamFileParser:

    def __init__(self, param_input) -> None:
        self.param_input = param_input
    
    def parse(self):
        param_dict = {}
        with open(self.param_input, 'r', encoding='utf-8') as file:
            lines = file.readlines()
    
            for idx, line in enumerate(lines, 1):
                params = line.rstrip().split(' ')
                params = tuple(float(x) for x in params)
                param_dict[idx] = params
        return param_dict


class Spectrum:
    """Class to represent a spectrum"""

    def __init__(self, name, lib_id=None, mw=None, precursor_mz=None, status=None, full_name=None, comment=None, peaks: List = None):
        self.name = name
        self.lib_id = lib_id
        self.mw = mw
        self.precursor_mz = precursor_mz
        self.status = status
        self.full_name = full_name
        self.comment = comment
        self.peaks = peaks

    def __str__(self):
        return f"Spectrum: {self.name}, PrecursorMZ: {self.precursor_mz}, NumPeaks: {len(self.peaks)}"


class SpTXTParser:
    
    def __init__(self) -> None:
        self.line = ""


    def parse_sptxt(self, filename):
        spectra = []
        current_spectrum = None

        with open(filename, 'r', encoding='utf-8') as file:
            for self.line in file:
                if self.line.startswith('Name: '):
                    if current_spectrum:
                        spectra.append(current_spectrum)
                    current_spectrum = self.parse_spectrum_header()
                elif current_spectrum and self.line.startswith(('LibID', 'MW', 'PrecursorMZ', 'Status', 'Comment')):
                    self.parse_spectrum_attribute(current_spectrum)
                elif current_spectrum and self.line.startswith('NumPeaks'):
                    num_peaks = int(self.line.split(':')[1].strip())
                    self.parse_spectrum_peaks(file, current_spectrum, num_peaks)

        return spectra


    def parse_spectrum_attribute(self, current_spectrum):
        attribute_name, attribute_value = map(str.strip, self.line.split(':', 1))
        attribute_value = float(attribute_value) if _is_numeric(attribute_value) else attribute_value
        setattr(current_spectrum, attribute_name.lower(), attribute_value)


    def parse_spectrum_peaks(self, file, current_spectrum, num_peaks):
        for _ in range(num_peaks):
            line = next(file).strip()
            self.parse_peak(line, current_spectrum)


    def parse_spectrum_header(self):
        name_match = re.match(r'Name: (.+)', self.line)
        if name_match:
            name = name_match.group(1)
        else:
            name = None

        return Spectrum(name=name)


    @staticmethod
    def parse_peak(peak_line, spectrum):
        parts = peak_line.split('\t')
        if len(parts) == 3:
            mz, intensity, annotations = parts
            spectrum.peaks.append((float(mz), float(intensity), annotations))


class FastaParser:

    def __init__(self, input_path) -> None:
        self.input_path = input_path
    

    def parse(self):

        current_id = None
        protein_dict = {}

        with open(self.input_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]
                    protein_dict[current_id] = ''
                elif current_id is not None:
                    protein_dict[current_id] += line

        return protein_dict
    

@dataclass
class PercolatorResult:
    psm_id: str
    score: float
    q_value: float
    posterior_error_prob: float
    peptide: str
    protein_ids: List



class PercolatorParser:
    
    def __init__(self, input_path) -> None:
        self.input_path = input_path


    def parse(self):

        with open(self.input_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

            results = [line.strip().split('\t') for line in lines[1:]] # skip first line with column names
            output_results = [PercolatorResult(x[0], *[float(i) for i in x[1:4]], x[4], x[5:]) for x in results]

        return output_results

