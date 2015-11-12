from . import unittest, THIS
import mock

class TestProjectRename(unittest.TestCase):
    def test_existing_imports_still_work(self):
        from bio_pieces import amos as bp_amos
        from bio_bits import amos as bb_amos
        self.assertEqual(bp_amos, bb_amos)
