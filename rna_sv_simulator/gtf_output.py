### TODO: Maybe adding parallel processing for each chromosome ???
import os


from mavis import constants
from mavis.annotate import file_io, fusion, variant
from mavis.breakpoint import Breakpoint, BreakpointPair

GTF_EMPTY = '.'  # value used to denote an empty field


def gtf_attr_to_string(attr_dict):
    '''
    Given some list of gtf attributes. Format them into a semi
    colon delimited string
    '''
    attributes_list = []
    for key, value in attr_dict.items():
        attributes_list.append('{} "{}"'.format(key, value))
    return ';'.join(attributes_list)


def output_gtf_file(annotations, output_path):
    '''
    Args:
        annotations (dict or list of Gene by str): list of genes keyed by chromosome name
        output_path (str): path the output annotations file
    '''

    header = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    row = {'source': 'rna_sv_simulator', 'attribute': {}}

    with open(output_path, 'w') as fh:
        for chrom in annotations:
            row['seqname'] = chrom
            for gene in annotations[chrom]:
                strand = gene.strand
                row['attribute'].update({
                    'gene_id': gene.name,
                    'gene_name': ','.join(gene.aliases),
                    'gene_start': gene.start,
                    'gene_end': gene.end
                })

                for pre_transcript in gene.transcripts:
                    row['attribute'].update({
                        'transcript_id': pre_transcript.name,
                        'transcript_start': pre_transcript.start,
                        'transcript_end': pre_transcript.end
                    })
                    for spliced_transcript in pre_transcript.spliced_transcripts:
                        for i, exon in enumerate(spliced_transcript.exons):
                            row.update({
                                'start': exon.start,
                                'end': exon.end
                            })

                            start_pos = str(exon.start)
                            end_pos = str(exon.end)
                            row['attribute'].update({
                                'exon_number': i + 1,
                                'exon_id': exon.name
                            })

                            row_copy = {}
                            row_copy.update(row)
                            row_copy['attribute'] = gtf_attr_to_string(row_copy['attribute'])

                            fh.write('\t'.join([
                                str(row_copy.get(col, GTF_EMPTY)) for col in header
                            ]) + '\n')

