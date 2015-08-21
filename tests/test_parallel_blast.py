import tempfile
import os

import mock

from bio_pieces.compat import unittest
from bio_pieces.compat import StringIO
from bio_pieces import parallel_blast

#http://talby.rcs.manchester.ac.uk/~ri/_notes_sge/par_envs_and_integration.html
# nodename cpucount queue processorrange
#$PE_HOSTFILE
SGE_HOSTFILE = '''node1.localhost 1 queue UNDEFINED
node2.localhost 2 queue UNDEFINED
node3.localhost 3 queue UNDEFINED
'''

#$PBS_NODEFILE
PBS_MACHINEFILE = '''node1.localhost
node2.localhost
node2.localhost
node3.localhost
node3.localhost
node3.localhost
'''

hosts = [
    ('node1.localhost', 1),
    ('node2.localhost', 2),
    ('node3.localhost', 3),
]
sshlogins = []
for x in hosts:
    sshlogins += ['--sshlogin', '{1}/{0}'.format(*x)]

class TestParseHostfile(unittest.TestCase):
    def test_parses_pbs_hostfile(self):
        r = parallel_blast.parse_hostfile(
            StringIO(PBS_MACHINEFILE)
        )
        self.assertListEqual(hosts, r)

    def test_parses_sge_hostfile(self):
        r = parallel_blast.parse_hostfile(
            StringIO(SGE_HOSTFILE)
        )
        self.assertListEqual(hosts, r)

    def test_unknown_hostfile(self):
        self.assertRaises(
            ValueError,
            parallel_blast.parse_hostfile,
            StringIO('node1.localhost something something something')
        )

class TestGetHostfile(unittest.TestCase):
    @mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PBS_NODEFILE': 'foo'})
    def test_detects_pbs(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('foo', r)
    
    @mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PE_HOSTFILE': 'foo'})
    def test_detects_sge(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('foo', r)

    def test_detects_local(self):
        r = parallel_blast.get_hostfile()
        self.assertEqual('', r)

class TestGenerateSSHLogins(unittest.TestCase):
    def setUp(self):
        _, self.hostfile = tempfile.mkstemp()
        self.addCleanup(os.unlink, self.hostfile)

    def test_single_compute_node_multiple_cpu(self):
        with open(self.hostfile, 'w') as fh:
            fh.write('node1.localhost\n')
            fh.write('node1.localhost\n')
        with mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(['--sshlogin', '2/:'], r)
        
    def test_pbs_sshlogins(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(PBS_MACHINEFILE)
        with mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_sge_sshlogins(self):
        with open(self.hostfile, 'w') as fh:
            fh.write(SGE_HOSTFILE)
        with mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PE_HOSTFILE': self.hostfile}):
            r = parallel_blast.generate_sshlogins()
            self.assertListEqual(sshlogins, r)

    def test_local_sshlogins(self):
        r = parallel_blast.generate_sshlogins(2)
        self.assertListEqual(['--sshlogin', '2/:'], r)

    def test_local_sshlogins_noninst(self):
        r = parallel_blast.generate_sshlogins()
        self.assertListEqual(['--sshlogin', ':'], r)

class MockSH(unittest.TestCase):
    def setUp(self):
        _, self.hostfile = tempfile.mkstemp()
        self.patch_sh_cmd = mock.patch('bio_pieces.parallel_blast.sh.Command')
        self.patch_sh_which = mock.patch('bio_pieces.parallel_blast.sh.which')
        self.mock_sh_which = self.patch_sh_which.start()
        self.mock_sh_cmd = self.patch_sh_cmd.start()
        self.addCleanup(self.patch_sh_cmd.stop)
        self.addCleanup(self.patch_sh_which.stop)
        self.patch_open = mock.patch('bio_pieces.parallel_blast.open')
        self.mock_open = self.patch_open.start()
        self.addCleanup(self.patch_open.stop)
        self.infile = '/path/infile'
        self.outfile = '/path/outfile'

class TestParallelBlast(MockSH):
    def test_correct_input_file_handling(self):
        self.mock_sh_which.return_value = '/path/to/foon'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        r = self.mock_sh_cmd.return_value.call_args
        # It seems that parallel needs 
        self.assertEqual(r[1]['_in'], self.mock_open.return_value)
        self.assertEqual(r[1]['_out'], self.mock_open.return_value)
        self.assertIn('--pipe', r[0])

    def test_handles_blastn_task(self):
        self.mock_sh_which.return_value = '/path/to/blastn'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'blastn', 'megablast',
            '-evalue 0.01 -otherblast arg'
        )
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('/path/to/blastn', r)
        self.assertIn('-task', r)
        self.assertIn('megablast', r)

    def test_handles_diamond(self):
        self.mock_sh_which.return_value = '/path/to/diamond'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'diamond', None,
            '-evalue 0.01 -otherblast arg'
        )
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('/path/to/diamond', r)
        self.assertNotIn('-task', r)

    def test_handles_blastx(self):
        self.mock_sh_which.return_value = '/path/to/blastx'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'blastx', None,
            '-evalue 0.01 -otherblast arg'
        )
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('/path/to/blastx', r)
        self.assertNotIn('-task', r)

    def test_command_string_is_correct(self):
        self.mock_sh_which.return_value = '/path/to/foon'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        self.mock_sh_cmd.assert_called_once_with('parallel')
        r = self.mock_sh_cmd.return_value.call_args
        blastcmd = r[0]
        print(r[0])
        self.assertIn('-db', blastcmd)
        self.assertIn('/path/db/nt', blastcmd)
        self.assertIn('-otherblast', blastcmd)
        self.assertIn('arg', blastcmd)
        self.assertIn('-max_target_seqs', blastcmd)
        self.assertIn(str(parallel_blast.MAX_TARGET_SEQS), blastcmd)
        self.assertIn('-outfmt', blastcmd)
        self.assertIn('"{0}"'.format(parallel_blast.BLAST_FORMAT), blastcmd)

    def test_localhost(self):
        self.mock_sh_which.return_value = '/path/to/foon'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        self.mock_sh_cmd.assert_called_once_with('parallel')
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('--sshlogin', r)
        self.assertIn('5/:', r)

    def test_remote_hosts(self):
        self.mock_open.return_value.__enter__.return_value = PBS_MACHINEFILE.splitlines()
        with mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            self.mock_sh_which.return_value = '/path/to/foon'
            parallel_blast.parallel_blast(
                self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
                '-evalue 0.01 -otherblast arg'
            )
            self.mock_sh_cmd.assert_called_once_with('parallel')
            r = self.mock_sh_cmd.return_value.call_args[0]
            self.assertEqual(3, r.count('--sshlogin'))
            self.assertIn('1/node1.localhost', r)
            self.assertIn('2/node2.localhost', r)
            self.assertIn('3/node3.localhost', r)

    def test_cannot_find_exe_raises_exception(self):
        self.mock_sh_which.return_value = None
        self.assertRaises(
            ValueError,
            parallel_blast.parallel_blast,
            self.infile, self.outfile, 5, '/path/to/db', 'blastn', 'foox', '-bar foo'
        )

class TestParallelDiamond(MockSH):
    def test_correct_inputoutput_handling(self):
        self.mock_sh_which.return_value = '/path/to/diamond'
        parallel_blast.parallel_blast(
            self.infile, self.outfile, 5, '/path/db/nt', 'foon', 'barn',
            '-evalue 0.01 -otherblast arg'
        )
        r = self.mock_sh_cmd.return_value.call_args
        # It seems that parallel needs 
        self.assertEqual(r[1]['_in'], self.mock_open.return_value)
        self.assertEqual(r[1]['_out'], self.mock_open.return_value)
        self.assertIn('--pipe', r[0])

    def test_cannot_find_exe_raises_exception(self):
        self.mock_sh_which.return_value = None
        self.assertRaises(
            ValueError,
            parallel_blast.parallel_diamond,
            self.infile, self.outfile, 5, '/path/to/dmd', 'foox', '-bar foo'
        )

    def test_correct_command_string(self):
        self.mock_sh_which.return_value = '/path/to/diamond'
        parallel_blast.parallel_diamond(
            self.infile, self.outfile, 5, '/path/to/dmd', 'foox', '-bar foo'
        )
        r = self.mock_sh_cmd.return_value.call_args[0]
        self.assertIn('foox', r)
        self.assertIn('--threads', r)
        self.assertIn('5', r)
        self.assertIn('--db', r)
        self.assertIn('/path/to/dmd', r)
        self.assertIn('--query', r)
        self.assertIn('{}', r)
        self.assertIn('--cat', r)
        self.assertIn('--sshlogin', r)
        self.assertIn('1/:', r)

    def test_each_remote_host_has_one_instance(self):
        self.mock_sh_which.return_value = '/path/to/diamond'
        self.mock_open.return_value.__enter__.return_value = PBS_MACHINEFILE.splitlines()
        with mock.patch.dict('bio_pieces.parallel_blast.os.environ', {'PBS_NODEFILE': self.hostfile}):
            parallel_blast.parallel_diamond(
                self.infile, self.outfile, 5, '/path/to/dmd', 'foox', '-bar foo'
            )
            r = self.mock_sh_cmd.return_value.call_args[0]
            self.assertIn('--sshlogin', r)
            self.assertIn('1/node1.localhost', r)
            self.assertIn('1/node2.localhost', r)
            self.assertIn('1/node3.localhost', r)

@mock.patch('bio_pieces.parallel_blast.parallel_blast')
@mock.patch('bio_pieces.parallel_blast.parallel_diamond')
class TestMain(unittest.TestCase):
    def setUp(self):
        self.patch_args = mock.patch('bio_pieces.parallel_blast.argparse.ArgumentParser')
        self.mock_args = self.patch_args.start()
        self.addCleanup(self.patch_args.stop)
        self.mock_args = self.mock_args.return_value.parse_args.return_value
        self.mock_args.outfile = '/path/out.blast'
        self.mock_args.ninst = 1
        self.mock_args.db = '/path/db'
        self.mock_args.task = None
        self.mock_args.blast_options = ''

    @mock.patch('bio_pieces.parallel_blast.exists', mock.Mock(return_value=True))
    def test_runs_parallel_diamond_for_diamond(self, mock_pdmnd, mock_pblast):
        self.mock_args.inputfasta = '/path/in.fa'
        self.mock_args.blast_exe = 'diamond'
        parallel_blast.main()
        self.assertEqual(1, mock_pdmnd.call_count)

    @mock.patch('bio_pieces.parallel_blast.exists', mock.Mock(return_value=True))
    def test_runs_parallel_blast_for_blast(self, mock_pdmnd, mock_pblast):
        self.mock_args.inputfasta = '/path/in.fa'
        self.mock_args.blast_exe = 'blastn'
        parallel_blast.main()
        self.assertEqual(1, mock_pblast.call_count)

    @mock.patch('bio_pieces.parallel_blast.exists', mock.Mock(return_value=False))
    def test_raises_assertion_when_input_fasta_missing(self, mock_pdmnd, mock_pblast):
        self.assertRaises(
            AssertionError,
            parallel_blast.main
        )
