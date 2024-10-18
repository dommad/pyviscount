import unittest
import configparser
import pandas as pd
import numpy as np
from partipy.postsearch import PostSearchValidation, PostSearchPartition
from partipy.template import TdcFdpFdrCalculator

class TestPostSearchPartition(unittest.TestCase):

    def setUp(self):
        # Create a sample DataFrame for testing
        data = {
            'sequence': ['AAA', 'BBB', 'CCC', 'AAA', 'BBB', 'CCC'],
            'modifications': ['mod1', 'mod2', 'mod3', 'mod1', 'mod2', 'mod3']
        }
        self.target_df = pd.DataFrame(data)


    def test_add_peptide_modification_index_to_target_df(self):
        # Create a sample instance of DataFrameProcessing

        # Call the add_peptide_modification_index_to_target_df method
        result_df, pep_mod_dict = PostSearchPartition(2).add_peptide_modification_index_to_target_df(self.target_df)

        # Perform assertions based on your expectations
        expected_columns = ['sequence', 'modifications', 'pep_mod', 'pep_mod_index']
        self.assertCountEqual(result_df.columns, expected_columns)

        # Replace with the expected pep_mod_dict based on your test case
        expected_pep_mod_dict = {('AAA', 'mod1'): 0, ('BBB', 'mod2'): 1, ('CCC', 'mod3'): 2}
        self.assertDictEqual(pep_mod_dict, expected_pep_mod_dict)


    def test_create_peptide_subset_id_mapping(self):
        # Create a sample instance of PostSearchPartition

        self.target_df['pep_mod_index'] = [0, 1, 2, 0, 1, 2]
        post_search_partition = PostSearchPartition(num_partitions=2)

        # Call the create_peptide_subset_id_mapping method
        result_mapping = post_search_partition.create_peptide_subset_id_mapping(self.target_df)

        # Perform assertions based on your expectations
        self.assertIsInstance(result_mapping, np.ndarray)
        self.assertEqual(len(result_mapping), len(self.target_df['pep_mod_index'].unique()))


class TestPostSearchValidation(unittest.TestCase):

    def setUp(self):
        # Create a sample DataFrame for testing
        data = pd.DataFrame({
            'scan': [1, 1, 2, 2, 3],
            'hit_rank': [1, 2, 1, 2, 1],
            'sequence': ['AAA', 'BBB', 'AAA', 'CCC', 'DDD'],
            'modifications': ['MOD1', 'MOD2', 'MOD1', 'MOD3', 'MOD4'],
            'tev': [0.18, 0.15, 0.15, 0.13, 0.33]
        })

        self.peptide_id_mapping = np.array([1, 0, 0, 0])
        
        self.target_df = pd.DataFrame(data)
        config_file_path = "./config.ini"
        with open(config_file_path, 'r', encoding='utf-8') as config_file:
            self.config = configparser.ConfigParser()
            self.config.read_file(config_file)

    

    def test_remove_below_threshold(self):
        # Create a sample instance of PostSearchValidation
        fdp_fdr_calculator = TdcFdpFdrCalculator  # Replace with actual initialization
        pep_mod_dict = {('AAA', 'MOD1'): 0, ('BBB', 'MOD2'): 1, ('CCC', 'MOD3'): 2, ('DDD', 'MOD4'): 3}  # Sample data, replace with actual values
        
        self.target_df['pep_mod_index'] = [0, 1, 0, 2, 3]
        post_search_validation = PostSearchValidation(self.config, fdp_fdr_calculator=fdp_fdr_calculator,
                                                      peptide_id_mapping=self.peptide_id_mapping, pep_mod_dict=pep_mod_dict,
                                                      target_df=self.target_df, td_df=pd.DataFrame(), decoy_df=pd.DataFrame())
        
        post_search_validation.subset_id = 0

        # Call the _remove_below_threshold method
        threshold = 0.17  # Replace with the actual threshold value
        result_indices = post_search_validation._remove_below_threshold(self.target_df, threshold)
        # Perform assertions based on your expectations
        
        expected_indices = np.array([1, 3])  # Replace with the expected indices based on your test case

        self.assertEqual(True, np.all(expected_indices == result_indices.values))



# class TestReadSearchResults(unittest.TestCase):
#     def setUp(self):
#         # Set up your test data or mocks here if needed
#         # Replace placeholders with your actual file paths and extensions
#         self.target_file = '../examples/adult_f16_short.pep.xml'
#         self.td_file = '../examples/adult_f16_short.pep.xml'
#         self.decoy_file = '../examples/adult_f16_short.pep.xml'
#         self.engine = 'Comet'  # Replace with your actual engine
#         self.threshold_score = 'tev'
#         self.fdr_score = 'tev'
#         self.test_instance = ps.PostSearchPartition(self.engine, self.threshold_score, self.fdr_score)

#         self.test_instance.target_df = pd.DataFrame({
#             'scan': [1, 1, 2, 2, 3],
#             'sequence': ['AAA', 'BBB', 'AAA', 'CCC', 'DDD'],
#             'modifications': ['MOD1', 'MOD2', 'MOD1', 'MOD3', 'MOD4'],
#             'tev': [0.7, 0.5, 0.8, 0.6, 0.9]
#         })



#     def test_read_search_results(self):
#         # Replace placeholders with your actual file paths and extensions
        
#         result = self.test_instance._read_search_results(self.target_file, self.td_file, self.decoy_file)

#         # Check if the result is a list of DataFrames
#         self.assertIsInstance(result, list)
#         for df in result:
#             self.assertIsInstance(df, pd.DataFrame)



#     def test_calculate_tev_score(self):
#         # Replace placeholders with your actual data or mocks
#         test_data = pd.DataFrame({'e_value': [1.0, 2.0, 3.0], 'p_value': [0.1, 0.2, 0.3], 'num_candidates': [10, 20, 30]})

#         # Replace placeholders with actual parameter values
#         parameter_a = 0.5
#         parameter_n0 = 100.0

#         result = self.test_instance.calculate_tev(test_data, parameter_a, parameter_n0)

#         # Check if the result is a pandas series
#         self.assertIsInstance(result, pd.core.series.Series)



#     def test_generate_peptide_id_mapping(self):
#         # Create an instance of SearchResultsProcessor
#         #test_instance = ps.PostSearchPartition(self.target_file, self.td_file, self.decoy_file, self.engine, self.file_ext)

        

#         mock_peptide_id_mapping = {('AAA', 'MOD1'): 0, ('BBB', 'MOD2'): 1, ('CCC', 'MOD3'): 1, ('DDD', 'MOD4'): 0}


#         # Call the method being tested
#         self.test_instance._create_peptide_subset_id_mapping()
#         result = self.test_instance.peptide_id_mapping

#         # Assert the type and structure of the result
#         self.assertIsInstance(result, dict)

#         # Add more specific assertions based on your expected output
#         # Ensure the keys are tuples and values are integers
#         for key, value in result.items():
#             self.assertIsInstance(key, tuple)
#             self.assertIsInstance(value, np.int64)
#             self.assertIn(key, mock_peptide_id_mapping.keys())



#     def test_remove_below_threshold(self):
#         # Create an instance of SearchResultsProcessor
#         #test_instance = ps.PostSearchPartition(self.target_file, self.td_file, self.decoy_file, self.engine, self.file_ext)
        

#         self.test_instance.peptide_id_mapping = {('AAA', 'MOD1'): 0, ('BBB', 'MOD2'): 1, ('CCC', 'MOD3'): 1, ('DDD', 'MOD4'): 0}

#         # Call the method being tested
#         result = self.test_instance._remove_below_threshold(threshold=0.75)

#         # Assert the type and structure of the result
#         self.assertIsInstance(result, list)

#         # Add more specific assertions based on your expected output
#         # Ensure the result contains the correct indices
#         self.assertListEqual(result, [2, 4])



#     def test_filter_update_td_df(self):
#         # Create an instance of SearchResultsProcessor
#         #test_instance = ps.PostSearch(self.target_file, self.td_file, self.decoy_file, self.engine, self.file_ext)
        

#         # Mock data for testing
#         mock_target_df = pd.DataFrame({
#             'scan': [1, 2, 3],
#             'sequence': ['AAA', 'BBB', 'CCC'],
#             'modifications': ['MOD1', 'MOD2', 'MOD3'],
#         })

#         self.test_instance.td_df = pd.DataFrame({
#             'scan': [1, 2, 3, 4],
#             'sequence': ['AAA', 'BBB', 'CCC', 'DDD'],
#             'modifications': ['MOD1', 'MOD2', 'MOD3', 'MOD4'],
#             'is_decoy': [False, False, False, True],
#             'hit_rank': [1, 2, 1, 1],
#         })

#         self.test_instance.peptide_id_mapping = {('AAA', 'MOD1'): 0, ('BBB', 'MOD2'): 1, ('CCC', 'MOD3'): 0}

#         # Call the method being tested
#         result = self.test_instance._filter_update_target_decoy_df(mock_target_df)

#         # Assert the type and structure of the result
#         self.assertIsInstance(result, pd.DataFrame)

#         # Add more specific assertions based on your expected output
#         # Ensure the DataFrame has the correct columns and values
#         self.assertListEqual(result.columns.tolist(), ['scan', 'sequence', 'modifications', 'is_decoy', 'hit_rank', 'gt_status'])
#         self.assertEqual(result['gt_status'].sum(), 2)


    
#     def test_calculate_fdp_fdr_contour(self):

#         self.test_instance.peptide_id_mapping = {('AAA', 'MOD1'): 0, ('BBB', 'MOD2'): 1, ('CCC', 'MOD3'): 1, ('DDD', 'MOD4'): 0}


#         self.test_instance.td_df = pd.DataFrame({
#             'scan': [1, 2, 3],
#             'tev': [0.3, 0.4, 0.2],
#             'sequence': ['AAA', 'BBB', 'CCC'],
#             'modifications': ['MOD1', 'MOD2', 'MOD3'],
#             'is_decoy': [False, False, False],
#             'hit_rank': [1, 2, 1],
#             'gt_status': [True, False, True]
#         })

#         # Call the method being tested
#         result = self.test_instance.calculate_fdp_fdr_contour()

#         # Assert the type and structure of the result
#         self.assertIsInstance(result, tuple)
#         self.assertEqual(len(result), 2)
#         self.assertIsInstance(result[0], list)
#         self.assertIsInstance(result[1], np.ndarray)

#         # Add more specific assertions based on your expected output
#         # ...



#     def test_calculate_proxy_fdp_fdr(self):
        
#         # Mock data for testing
#         mock_df = pd.DataFrame({
#             'is_decoy': [False, False, True, True, False],
#             'gt_status': [True, False, False, False, True],
#             'tev': [0.8, 0.6, 0.7, 0.9, 0.5],
#         })

#         # Call the method being tested
#         result = self.test_instance.calculate_proxy_fdp_fdr(mock_df)

#         # Assert the type and structure of the result
#         self.assertIsInstance(result, tuple)
#         self.assertEqual(len(result), 2)
#         self.assertIsInstance(result[0], np.ndarray)
#         self.assertIsInstance(result[1], np.ndarray)

#         # Add more specific assertions based on your expected output
#         # ...




if __name__ == '__main__':
    unittest.main()