import random
import unittest
from unittest.mock import MagicMock, patch

from mavis.annotate import genomic as _genomic
from mavis import breakpoint as _breakpoint
from mavis import constants

from rna_sv_simulator import fusion

class TestMutateContinuous(unittest.TestCase):

    def setUp(self):
        self.chr1 = ''.join([random.choice('ATCG') for i in range(1000)])

        self.reference_genome = {
            '1': MagicMock(seq=self.chr1),
        }

        self.annotations = {'1': [
            self.create_gene('1', '+', [(100, 200), (250, 300), (320, 360)]),
            self.create_gene('1', '-', [(800, 820), (830, 850)]),
            self.create_gene('1', '+', [(300, 350), (420, 500)]),
            self.create_gene('1', '-', [(300, 350), (380, 500)]),
            self.create_gene('1', '-', [(800, 820), (830, 850)]),
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
    
    def mutate(self, strand1, strand2, orient1, orient2, event):
        bpp = _breakpoint.BreakpointPair(
            _breakpoint.Breakpoint('1', 400, strand=strand1, orient=orient1),
            _breakpoint.Breakpoint('1', 700, strand=strand2, orient=orient2),
            untemplated_seq='',
            event_type=event
        )

        new_reference_genome, new_annotations = fusion._mutate_continuous(
            self.reference_genome,
            self.annotations,
            bpp)
        return bpp, new_reference_genome, new_annotations

    # deletion event type
    def test_del_lr_pos_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '+', 'L', 'R', constants.SVTYPE.DEL
        )
        exp_seq = self.chr1[:400] + self.chr1[699:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_del_lr_neg_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '-', 'L', 'R', constants.SVTYPE.DEL
        )
        exp_seq = self.chr1[:400] + self.chr1[699:]
        self.assertEqual(len(exp_seq), len(new_reference_genome))
        self.assertEqual(exp_seq, new_reference_genome)


    # duplication event type
    def test_dup_rl_pos_pos(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '+', 'R', 'L', constants.SVTYPE.DUP
        )
        exp_seq = self.chr1[:399] + self.chr1[399:700] + self.chr1[399:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_dup_rl_neg_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '-', 'R', 'L', constants.SVTYPE.DUP
        )
        exp_seq = self.chr1[:399] + self.chr1[399:700] + self.chr1[399:]
        self.assertEqual(exp_seq, new_reference_genome)


    # inversion event type
    def test_inv_ll_pos_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '-', 'L', 'L', constants.SVTYPE.INV
        )
        exp_seq = self.chr1[:400] + constants.reverse_complement(self.chr1[400:700]) + self.chr1[700:]
        self.assertEqual(len(exp_seq), len(self.chr1))
        self.assertEqual(exp_seq, new_reference_genome)

    def test_inv_ll_neg_post(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '+', 'L', 'L', constants.SVTYPE.INV
        )
        exp_seq = self.chr1[:400] + constants.reverse_complement(self.chr1[400:700]) + self.chr1[700:]
        diffs = []
        for (i, (exp, act)) in enumerate(zip(exp_seq, new_reference_genome)):
            if (exp != act):
                diffs.append((i, (exp, act, self.chr1[i])))
        print(diffs)
        self.assertEqual(exp_seq, new_reference_genome)

    def test_inv_rr_pos_neg(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '+', '-', 'R', 'R', constants.SVTYPE.INV
        )
        exp_seq = self.chr1[:399] + constants.reverse_complement(self.chr1[399:699]) + self.chr1[699:]
        self.assertEqual(exp_seq, new_reference_genome)

    def test_inv_rr_neg_post(self):
        bpp, new_reference_genome, new_annotations = self.mutate(
            '-', '+', 'R', 'R', constants.SVTYPE.INV
        )
        exp_seq = self.chr1[:399] + constants.reverse_complement(self.chr1[399:699]) + self.chr1[699:]
        self.assertEqual(exp_seq, new_reference_genome)

    


