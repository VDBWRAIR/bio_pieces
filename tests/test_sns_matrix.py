from . import unittest

from bio_pieces import sns_matrix

class TestSNSMatrix(unittest.TestCase):
    def test_sns_matrix_same_nt_returns_10(self):
        for c in 'ACGINRTUXY':
            self.assertEqual(10, sns_matrix.sns_matrix[c][c])

    def test_sns_matrix_returns_correct_values(self):
        self.assertEqual(-8, sns_matrix.sns_matrix['A']['C'])
        self.assertEqual(-8, sns_matrix.sns_matrix['C']['A'])
