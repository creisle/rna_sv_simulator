'''
create the fusion chromosome and re-mapping of the genes on the mutant chrs
'''
from mavis.constants import SVTYPE, ORIENT, STRAND
from mavis.annotate import fusion as _fusion
from mavis.annotate import genomic as _genomic


def determine_breakpoint_prime(breakpoint_pair):
    '''
    For a given breakpoint return True if break1 is five prime
    for the mutant chr and return false if break2 is five prime
    for the mutant chr
    '''
    break1 = breakpoint_pair.break1
    break2 = breakpoint_pair.break2

    if breakpoint_pair.opposing_strands:
        # orientations are the same
        if break1.orient == ORIENT.LEFT:
            if break2.strand == STRAND.POS:
                return False
        elif break2.strand == STRAND.NEG:
            return False
    else:
        # strands are the same
        if break1.strand == STRAND.POS:
            if break2.orient == ORIENT.LEFT:
                return False
        elif break2.orient == ORIENT.RIGHT:
            return False
    return True


def get_chr_seq(breakpoint, reference_genome):
    seq = ''
    if breakpoint.orient == ORIENT.LEFT:
        seq = reference_genome[breakpoint.chr].seq[:breakpoint.end]
    else:
        seq = reference_genome[breakpoint.chr].seq[breakpoint.start:]

    if breakpoint.strand == STRAND.NEG:
        seq = reverse_complement(seq)
    return seq


def _mutate_translocation(reference_genome, annotations, breakpoint_pair):

    new_chrom_refseq = ''
    new_annotations = []

    break1 = breakpoint_pair.break1
    break2 = breakpoint_pair.break2

    break1_is_five_prime = determine_breakpoint_prime(breakpoint_pair)
    break1_offset = 0
    break2_offset = 0

    # create the new sequence
    if break1_is_five_prime:
        new_chrom_refseq = get_chr_seq(break1, reference_genome) + get_chr_seq(break2, reference_genome)
    else:
        new_chrom_refseq = get_chr_seq(break2, reference_genome) + get_chr_seq(break1, reference_genome)

    # shift all the gene annotations based on the new chr
    if not break1_is_five_prime:
        break1, break2 = break2, break1

    for gene in annotations[break1chr]:
        if break1.orient == ORIENT.LEFT:
            new_annotations.append(gene)
        else:
            new_annotations.append(
                shift_gene(gene, lambda x: len(reference_genome[break1.chr].seq) - x)
            )



def shift_gene(gene, offset_func, flipped=False):
    '''

    '''
    new_strand = gene.strand
    if flipped:
        new_strand = STRAND.POS if gene.strand == STRAND.NEG else STRAND.NEG
    new_gene = _genomic.Gene(
        gene.chr,
        offset_func(gene.end if flipped else gene.start),
        offset_func(gene.start if flipped else gene.end),
        name=gene.name,
        strand=new_strand,
        aliases=gene.aliases
    )
    for transcript in gene.unspliced_transcripts:
        new_exons = []
        for exon in transcript.exons:
            new_exon = _genomic.Exon(
                offset_func(gene.end if flipped else gene.start),
                offset_func(gene.start if flipped else gene.end),
                strand=new_strand,
                name=exon.name
            )
            new_exons.append(exon)
        new_pre_transcript = _genomic.PreTranscript(
            exons=new_exons,
            gene=new_gene,
            name=transcript.name,
            is_best_transcript=transcript.is_best_transcript
        )
        new_gene.unspliced_transcripts.append(new_pre_transcript)
        for splicing_pattern in new_pre_transcript.generate_splicing_patterns():
            new_pre_transcript.transcripts.append(
                _genomic.Transcript(
                    new_pre_transcript,
                    splicing_pattern
                )
            )
    return new_gene


def _mutate_continuous(reference_genome, annotations, breakpoint_pair):
    # create the original fusion transcirpt
    ft = _fusion.FusionTranscript.build()

    mutant_genes = []
    mutant_seq = ''

    for gene in annotations[breakpoint_pair.break1.chr]:
        if gene.end < breakpoint_pair.break1.start:
            # all genes to the left remain the same
            mutant_genes.append(gene)
        elif gene.start > breakpoint_pair.break2.end:
            # genes to the right are affected by the event
            mutant_gene = shift_gene(gene, ...) # TODO
        else:
            # genes in the middle depend on the event type
            pass


def mutate(reference_genome, annotations, breakpoint_pair):
    '''
    Returns the mutated reference genome sequences as well as the
    shifted annotations and the fusion transcript
    '''
    if breakpoint_pair.interchromosomal:
        return _mutate_translocation()
    else:
        return _mutate_continuous()


