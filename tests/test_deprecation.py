from . import unittest, THIS
import mock
import sys

dep_msg = 'bio_pieces will be removed in the future. Please use import bio_bits'

@mock.patch('bio_bits.version.get_release', mock.MagicMock(return_value=1.0))
@mock.patch.object(sys, 'stderr')
class TestProjectRename(unittest.TestCase):
    def test_can_import_bio_pieces(self, mock_serr):
        import bio_pieces
        r = bio_pieces.version.get_release()
        self.assertEqual(1.0, r)
        actual_msg = mock_serr.write.call_args[0][0]
        self.assertIn(dep_msg, actual_msg)

