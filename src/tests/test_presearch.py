import unittest
import pandas as pd
from partipy.presearch import PreSearchPartition
from partipy.template import TdcFdpFdrCalculator as tdc

class TestAddPeptideModificationIndex(unittest.TestCase):

    def setUp(self) -> None:

        data = {'sequence': ['ABC', 'DEF', 'GHI'],
                'modifications': ['mod1', 'mod2', 'mod1']}
        self.df = pd.DataFrame(data)

        self.presearch_instance = PreSearchPartition(self.df, self.df, self.df, 'tev', 1, 'tev', 5, tdc)


    def test_add_peptide_modification_index(self):
        # Create a sample DataFrame for testing

        # Call the function
        result_df, index_mapping = self.presearch_instance.add_peptide_modification_index_to_target_df(self.df)

        # Assert that the DataFrame has the expected columns
        self.assertTrue('pep_mod' in result_df.columns)
        self.assertTrue('pep_mod_idx' in result_df.columns)

        



if __name__ == '__main__':
    unittest.main()
