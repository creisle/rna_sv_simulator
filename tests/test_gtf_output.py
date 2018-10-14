import unittest
from unittest import mock, patch, mock_open
from rna_sv_simulator import gtf_output

from mavis.annotate import genomic as _genomic
from mavis import breakpoint as _breakpoint
from mavis import constants


class TestOutputGtfFile(unittest.TestCase):

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

    def test_gtf_output(self):
        m = mock_open()
        with patch("builtins.open", m):
            gtf_output.output_gtf_file(self.annotations, '../letters')
        m.assert_called_with('../letters', 'w')
        handle_file = m()
        handle_file.write.assert_called_with('12\tprocessed_transcript\texon\t175931\t176602\t.\t+\t.\tgene_id\t"ENSG00000120645";gene_name "IQSEC3";gene_start "175931";gene_end "287626";transcript_id "ENST00000538872";transcript_start "175931";transcript_end "287626";exon_index "1"')
