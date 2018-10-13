### Maybe adding parallel processing for each chromosome ???
from mavis import constants
from mavis.annotate import file_io, fusion, variant
from mavis.breakpoint import Breakpoint, BreakpointPair

def attribute_dict_to_str(attr_dict):
    attributes_list = []
    for key,values in attr_dict.items():
        attributes_list.append(key + ' ' + '"' + values + '"')

    return ";".join(attributes_list)

def output_gtf_file(annot_dict, output_path):
    chromosome_id = '' # keeps the chromosome ids
    trans_name = 'processed_transcript'
    feature_type = 'exon' # Static since we are going to include all exons
    start_pos = '' # Exon start position
    end_pos = '' # Exon end position
    score = '.'
    strand = ''  # Exon Strand coming from gene_strand
    frame = '.'
    attribute = {} # This variable will be consists other variables to keep in gtf
    file_ = open(output_path+"result.gtf", "w") # file name should be specified
    for chromosome in annot_dict.keys():
        chromosome_id = chromosome
        for gene in annot_file_json[chrom]:
            strand = gene.strand
            attribute['gene_id'] = gene.name
            attribute['gene_name'] = ",".join(gene.aliases)
            attribute['gene_start'] = str(gene.start)
            attribute['gene_end'] = str(gene.end)

            for transcript in gene.transcripts:
                attribute['transcript_id'] = transcript.name
                attribute['transcript_start'] = str(transcript.start)
                attribute['transcript_end'] = str(transcript.end)

                for spliced_transcript in transcript.spliced_transcripts:
                    itr_cnt = 1
                    for exon in spliced_transcript.exons:
                        start_pos = str(exon.start)
                        end_pos = str(exon.end)
                        attribute['exon_index'] = str(itr_cnt)
                        itr_cnt += 1

                        ## Writing it to a file:
                        attribute_str = attribute_dict_to_str(attribute)

                        row_info = (chromosome_id + '\t' + trans_name + '\t' + feature_type + '\t'
                                    + start_pos + '\t' + end_pos + '\t' + score + '\t' + strand + '\t'
                                    + frame + '\t' + attribute_str + '\n')

                        file_.write(row_info) # Write to a file

    file_.close()
