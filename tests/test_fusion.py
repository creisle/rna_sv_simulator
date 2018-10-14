import random
import unittest
from unittest.mock import MagicMock, patch

from mavis.annotate import genomic as _genomic
from mavis import breakpoint as _breakpoint
from mavis import constants

from rna_sv_simulator import fusion


class TestShiftGenes(unittest.TestCase):

    def setUp(self):
        self.original_gene = _genomic.Gene('1', 1000, 2000, strand='+')

    def test_positive_shift(self):
        offset = lambda x: x + 10
        shifted_gene = fusion.shift_gene(self.original_gene, offset)

        self.assertEqual(self.original_gene.start + 10, shifted_gene.start)
        self.assertEqual(self.original_gene.end + 10, shifted_gene.end)

    def test_negative_shift(self):
        offset = lambda x: x - 10
        shifted_gene = fusion.shift_gene(self.original_gene, offset)

        self.assertEqual(self.original_gene.start - 10, shifted_gene.start)
        self.assertEqual(self.original_gene.end - 10, shifted_gene.end)

    def test_flip(self):
        offset = lambda x: 3000 - x
        shifted_gene = fusion.shift_gene(self.original_gene, offset, True)

        self.assertEqual(3000 - self.original_gene.start, shifted_gene.end)
        self.assertEqual(3000 - self.original_gene.end, shifted_gene.start)


def mock_rc(s):
    return s[::-1]


class TestGetChrSeq(unittest.TestCase):

    def setUp(self):
        self.reference_genome = {
            '1': MagicMock(seq='abcdefghijklmnopqrstuvwxyz')
        }

    def test_left_pos(self):
        breakpoint = _breakpoint.Breakpoint('1', 6, strand='+', orient='L')
        seq = fusion.get_chr_seq(breakpoint, self.reference_genome)
        self.assertEqual('abcdef', seq)

    @patch.object(constants, 'reverse_complement', new=mock_rc)
    def test_left_neg(self):
        breakpoint = _breakpoint.Breakpoint('1', 6, strand='-', orient='L')
        seq = fusion.get_chr_seq(breakpoint, self.reference_genome)
        self.assertEqual('abcdef'[::-1], seq)

    def test_right_pos(self):
        breakpoint = _breakpoint.Breakpoint('1', 6, strand='+', orient='R')
        seq = fusion.get_chr_seq(breakpoint, self.reference_genome)
        self.assertEqual('fghijklmnopqrstuvwxyz', seq)

    @patch.object(constants, 'reverse_complement', new=mock_rc)
    def test_right_neg(self):
        breakpoint = _breakpoint.Breakpoint('1', 6, strand='-', orient='R')
        seq = fusion.get_chr_seq(breakpoint, self.reference_genome)
        self.assertEqual('fghijklmnopqrstuvwxyz'[::-1], seq)


class TestDetermineBreakpointPrime(unittest.TestCase):
    def test_rl_pos(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='+', orient='R'),
            _breakpoint.Breakpoint('1', 15, strand='+', orient='L')
        )
        self.assertFalse(fusion.determine_breakpoint_prime(bpp))

    def test_rl_neg(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='-', orient='R'),
            _breakpoint.Breakpoint('1', 15, strand='-', orient='L')
        )
        self.assertTrue(fusion.determine_breakpoint_prime(bpp))

    def test_lr_pos(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='+', orient='L'),
            _breakpoint.Breakpoint('1', 15, strand='+', orient='R')
        )
        self.assertTrue(fusion.determine_breakpoint_prime(bpp))

    def test_lr_neg(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='-', orient='L'),
            _breakpoint.Breakpoint('1', 15, strand='-', orient='R')
        )
        self.assertFalse(fusion.determine_breakpoint_prime(bpp))

    def test_ll_pos_neg(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='+', orient='L'),
            _breakpoint.Breakpoint('1', 15, strand='-', orient='L')
        )
        self.assertTrue(fusion.determine_breakpoint_prime(bpp))

    def test_rr_pos_neg(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='+', orient='R'),
            _breakpoint.Breakpoint('1', 15, strand='-', orient='R')
        )
        self.assertFalse(fusion.determine_breakpoint_prime(bpp))

    def test_ll_neg_pos(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='-', orient='L'),
            _breakpoint.Breakpoint('1', 15, strand='+', orient='L')
        )
        self.assertFalse(fusion.determine_breakpoint_prime(bpp))

    def test_rr_neg_pos(self):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 10, strand='-', orient='R'),
            _breakpoint.Breakpoint('1', 15, strand='+', orient='R')
        )
        self.assertTrue(fusion.determine_breakpoint_prime(bpp))


class TestMutateTranslocation(unittest.TestCase):

    def setUp(self):
        self.chr1 = ''.join([random.choice('ATCG') for i in range(1000)])
        self.chr2 = ''.join([random.choice('ATCG') for i in range(1000)])
        self.reference_genome = {
            '1': MagicMock(seq=self.chr1),
            '2': MagicMock(seq=self.chr2)
        }
        self.annotations = {'1': [
            self.create_gene('1', '+', [(100, 200), (250, 300), (380, 460)]),
            self.create_gene('1', '-', [(800, 820), (830, 850)]),
        ], '2': [
            self.create_gene('2', '-', [(270, 300), (310, 360), (420, 440)]),
            self.create_gene('2', '+', [(700, 720), (800, 830)]),
        ]}

    def create_gene(self, chr, strand, exons):
        start = min([x for x, y in exons]) - 10
        end = max([y for x, y in exons]) + 10
        gene = _genomic.Gene(chr, start, end, strand=strand)
        new_exons = [_genomic.Exon(start, end, strand=strand) for start, end in exons]
        pre_transcript = _genomic.PreTranscript(new_exons, gene=gene)
        gene.unspliced_transcripts.append(pre_transcript)
        for spl_patt in pre_transcript.generate_splicing_patterns():
            transcript = _genomic.Transcript(pre_transcript, spl_patt)
            pre_transcript.transcripts.append(transcript)
        return gene

    def mutate(self, strand1, strand2, orient1, orient2):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 500, strand=strand1, orient=orient1),
            _breakpoint.Breakpoint('2', 500, strand=strand2, orient=orient2),
            untemplated_seq='A'
        )

        new_reference_genome, new_annotations = fusion._mutate_translocation(
            self.reference_genome,
            self.annotations,
            bpp)
        return bpp, new_reference_genome, new_annotations

    def test_itrans_ll_pos_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '-', 'L', 'L'
        )
        exp_seq = self.chr1[:500] + 'A' + constants.reverse_complement(self.chr2[:500])
        self.assertEqual(exp_seq, new_reference_genome)

    def test_itrans_rr_pos_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '-', 'R', 'R'
        )
        exp_seq = constants.reverse_complement(self.chr2[499:]) + 'T' + self.chr1[499:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_trans_lr_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '+', 'L', 'R'
        )
        exp_seq = self.chr1[:500] + 'A' + self.chr2[499:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_trans_rl_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '+', 'R', 'L'
        )
        exp_seq = self.chr2[:500] + 'A' + self.chr1[499:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_itrans_ll_neg_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '+', 'L', 'L'
        )
        exp_seq = self.chr2[:500] + 'A' + constants.reverse_complement(self.chr1[:500])
        self.assertEqual(exp_seq, new_reference_genome)

    def test_itrans_rr_neg_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '+', 'L', 'L'
        )
        exp_seq = self.chr2[:500] + 'A' + constants.reverse_complement(self.chr1[:500])
        self.assertEqual(exp_seq, new_reference_genome)

    def test_trans_lr_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '-', 'L', 'R'
        )
        exp_seq = '{}T{}'.format(
            constants.reverse_complement(self.chr2[499:]),
            constants.reverse_complement(self.chr1[:500])
        )
        self.assertEqual(exp_seq, new_reference_genome)

    def test_trans_rl_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '-', 'R', 'L'
        )
        exp_seq = '{}T{}'.format(
            constants.reverse_complement(self.chr1[499:]),
            constants.reverse_complement(self.chr2[:500])
        )
        self.assertEqual(exp_seq, new_reference_genome)


class TestGetReciprocal(unittest.TestCase):

    def create_bpp(self, break1, strand1, orient1, break2, strand2, orient2):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', break1, strand=strand1, orient=orient1),
            _breakpoint.Breakpoint('2', break2, strand=strand2, orient=orient2),
            untemplated_seq='A'
        )
        return bpp

    def test_itrans_ll_pos_neg(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '+', 'L', 500, '-', 'L'))
        exp = self.create_bpp(501, '+', 'R', 501, '-', 'R')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_itrans_rr_pos_neg(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '+', 'R', 500, '-', 'R'))
        exp = self.create_bpp(499, '+', 'L', 499, '-', 'L')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_trans_lr_pos(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '+', 'L', 500, '+', 'R'))
        exp = self.create_bpp(501, '+', 'R', 499, '+', 'L')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_trans_rl_pos(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '+', 'R', 500, '+', 'L'))
        exp = self.create_bpp(499, '+', 'L', 501, '+', 'R')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_itrans_ll_neg_pos(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '-', 'L', 500, '+', 'L'))
        exp = self.create_bpp(501, '-', 'R', 501, '+', 'R')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_itrans_rr_neg_pos(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '-', 'R', 500, '+', 'R'))
        exp = self.create_bpp(499, '-', 'L', 499, '+', 'L')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_trans_lr_neg(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '-', 'L', 500, '-', 'R'))
        exp = self.create_bpp(501, '-', 'R', 499, '-', 'L')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

    def test_trans_rl_neg(self):
        rep = fusion.get_reciprocal(self.create_bpp(500, '-', 'R', 500, '-', 'L'))
        exp = self.create_bpp(499, '-', 'L', 501, '-', 'R')
        self.assertEqual(exp.break1, rep.break1)
        self.assertEqual(exp.break2, rep.break2)

