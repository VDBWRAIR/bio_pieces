from . import (
    unittest, THIS
)

import shutil
import subprocess
import os
from os.path import *

import mock

from bio_bits import rename_fasta

TEST_FASTA = join(THIS,'example.fasta')
TEST_CSV = join(THIS,'renamelist.csv')

FASTA = '''>id1
ATGC
>id2
ATGC
'''

class TestGetCSVMap(unittest.TestCase):
    def test_gets_mapping(self):
        mapping = rename_fasta.get_csv_map(TEST_CSV)
        self.assertEqual(mapping['001'], 'sample1')
        self.assertEqual(mapping['002'], 'sample2')
        #self.assertEqual(mapping['003'], 'sample3')
        self.assertNotIn('#From', mapping)
        self.assertNotIn('003', mapping)

class TestRenameFastaIdentifiers(unittest.TestCase):
    def setUp(self):
        self.lines = FASTA.splitlines()
        self.patch_fileinput = mock.patch('bio_bits.rename_fasta.fileinput')
        self.mock_fileinput = self.patch_fileinput.start()
        self.mock_fileinput_input = self.mock_fileinput.input
        self.addCleanup(self.mock_fileinput.stop)
        self.patch_sys_stdout = mock.patch('bio_bits.rename_fasta.sys.stdout')
        self.patch_sys_stderr = mock.patch('bio_bits.rename_fasta.sys.stderr')
        self.mock_stdout = self.patch_sys_stdout.start()
        self.mock_stderr = self.patch_sys_stderr.start()
        self.addCleanup(self.mock_stdout.stop)
        self.addCleanup(self.mock_stderr.stop)

    def test_renames_to_console(self):
        self.mock_fileinput_input.return_value = self.lines
        mapping = {
            'id1': 'foo',
            'id2': 'bar'
        }
        rename_fasta.rename_fasta_identifiers(['foo.fasta'], mapping, False)
        self.mock_fileinput_input.assert_called_once_with(
            ['foo.fasta'], inplace=False
        )
        self.assertEqual(
            [
                mock.call('>foo'), mock.call('ATGC'),
                mock.call('>bar'), mock.call('ATGC')
            ],
            self.mock_stdout.write.call_args_list
        )

    def test_renames_inplace(self):
        self.mock_fileinput_input.__iter__ = []
        mapping = {
            'id1': 'foo',
            'id2': 'bar'
        }
        rename_fasta.rename_fasta_identifiers(['foo.fasta'], mapping, True)
        self.mock_fileinput_input.assert_called_once_with(
            ['foo.fasta'], inplace=True
        )

    def test_missing_identifier_in_mapping_file_stderr_message(self):
        self.mock_fileinput_input.return_value = self.lines
        mapping = {
            'id1': 'foo'
        }
        rename_fasta.rename_fasta_identifiers(['foo.fasta'], mapping, False)
        self.mock_fileinput_input.assert_called_once_with(
            ['foo.fasta'], inplace=False
        )
        self.assertEqual(
            [
                mock.call('>foo'), mock.call('ATGC'),
                mock.call('>id2'), mock.call('ATGC'),
            ],
            self.mock_stdout.write.call_args_list
        )
        self.mock_stderr.write.assert_called_once_with(
            'id2 is not in provided mapping\n'
        )

class TestFunctional(unittest.TestCase):
    def setUp(self):
        self.test_fasta = '/tmp/test.fasta'
        shutil.copyfile(TEST_FASTA, self.test_fasta)
        self.addCleanup(os.unlink, self.test_fasta)
        self.mapping = rename_fasta.get_csv_map(TEST_CSV)
        self.rev_mapping = {}
        for key,value in self.mapping.items():
            self.rev_mapping[value] = key
        self.rev_mapping['003'] = '003'

        self.patcher_argparse = mock.patch('bio_bits.rename_fasta.argparse')
        self.mock_argparse = self.patcher_argparse.start()
        self.addCleanup(self.mock_argparse.stop)

    def test_renames_inplace(self):
        cmd = [
            'rename_fasta',
            TEST_CSV,
            self.test_fasta,
            '--inplace'
        ]

        self._run_check_script(cmd, self.rev_mapping)

    def test_renames_stdin_stdout(self):
        cmd = [
            'rename_fasta',
            TEST_CSV,
            '-',
        ]

        sin = open(self.test_fasta)

        self._run_check_script(cmd, self.rev_mapping, sin)

    def _run_check_script(self, cmd, rev_mapping, stdin=None):
        if stdin is None:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                stdin=stdin
            )
        sout, serr = p.communicate()

        if '--inplace' in cmd:
            _iter = open(self.test_fasta)
        else:
            _iter = sout.splitlines()

        for line in _iter:
            line = str(line)
            if line.startswith('>'):
                p = line.split()
                id = p[0][1:]
                self.assertIn(id, rev_mapping)

        self.assertEqual(b'003 is not in provided mapping\n', serr)
