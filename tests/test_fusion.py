import unittest
from unittest.mock import MagicMock

from mavis.annotate import genomic as _genomic

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


class TestGetChrSeq(unittest.TestCase):

    def setUp(self):
        self.reference_genome = {
            '1': MagicMock(seq='abcdefghijklmnopqrstuvwxyz')
        }



