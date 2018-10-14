'''
create the fusion chromosome and re-mapping of the genes on the mutant chrs
'''
from mavis.constants import SVTYPE, ORIENT, STRAND
from mavis import constants
from mavis.annotate import fusion as _fusion
from mavis.annotate import genomic as _genomic
from mavis import breakpoint as _breakpoint


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
        seq = reference_genome[breakpoint.chr].seq[breakpoint.start - 1:]

    if breakpoint.strand == STRAND.NEG:
        seq = constants.reverse_complement(seq)
    return seq


def get_reciprocal(breakpoint_pair):
    '''
    Given a breakpoint pair representing a translocation
    return a breakpoint pair representing its reciprocal
    '''
    shift1 = 1 if breakpoint_pair.break1.orient == ORIENT.LEFT else -1
    break1 = _breakpoint.Breakpoint(
        breakpoint_pair.break1.chr,
        breakpoint_pair.break1.start + shift1,
        breakpoint_pair.break1.end + shift1,
        strand=breakpoint_pair.break1.strand,
        orient=ORIENT.LEFT if breakpoint_pair.break1.orient == ORIENT.RIGHT else ORIENT.RIGHT
    )
    shift2 = 1 if breakpoint_pair.break2.orient == ORIENT.LEFT else -1
    break2 = _breakpoint.Breakpoint(
        breakpoint_pair.break2.chr,
        breakpoint_pair.break2.start + shift2,
        breakpoint_pair.break2.end + shift2,
        strand=breakpoint_pair.break2.strand,
        orient=ORIENT.LEFT if breakpoint_pair.break2.orient == ORIENT.RIGHT else ORIENT.RIGHT
    )
    return _breakpoint.BreakpointPair(
        break1, break2,
        stranded=breakpoint_pair.stranded,
        untemplated_seq=breakpoint_pair.untemplated_seq,
        opposing_strands=breakpoint_pair.opposing_strands,
        **breakpoint_pair.data
    )


def _mutate_translocation(reference_genome, annotations, breakpoint_pair):
    new_annotations = []

    break1 = breakpoint_pair.break1
    break2 = breakpoint_pair.break2

    break1_is_five_prime = determine_breakpoint_prime(breakpoint_pair)

    # create the new sequence
    break1_seq = get_chr_seq(break1, reference_genome)
    break2_seq = get_chr_seq(break2, reference_genome)

    # shift all the gene annotations based on the new chr
    if not break1_is_five_prime:
        break1, break2 = break2, break1
        break1_seq, break2_seq = break2_seq, break1_seq

    new_annotations = []

    for gene in annotations[break1.chr]:
        # only move the genes retained and not covering the breakpoint
        if break1.orient == ORIENT.LEFT:   # positive strand
            if gene.end <= break1.start:
                new_gene = shift_gene(gene, lambda pos: pos)
                break1_seq += breakpoint_pair.untemplated_seq
            else:
                continue
        elif gene.start >= break1.end:   # negative strand
            offset = lambda pos: len(reference_genome[break1.chr].seq) - pos
            new_gene = shift_gene(gene, offset, True)
            break1_seq += constants.reverse_complement(breakpoint_pair.untemplated_seq)
        else:
            continue

        new_annotations.append(new_gene)

    for gene in annotations[break2.chr]:
        # only move the genes retained and not covering the breakpoint
        if break2.orient == ORIENT.RIGHT:  # positive strand
            if gene.start >= break2.end:
                offset = lambda pos: len(break1_seq) + pos - break2.end
                new_gene = shift_gene(gene, offset)
            else:
                continue
        elif gene.end <= break2.start:  # negative strand
            offset = lambda pos: len(break1_seq) + break2.end - pos
            new_gene = shift_gene(gene, offset, True)
        else:
            continue

        new_annotations.append(new_gene)

    return break1_seq + break2_seq, new_annotations


def shift_gene(gene, offset_func, flipped=False):
    '''
    Given some gene and a function calculating the shift in a genomic position
    return a new gene with positions of the gene and child exons moved
    '''
    if flipped:
        new_strand = STRAND.POS if gene.strand == STRAND.NEG else STRAND.NEG
    else:
        new_strand = gene.strand

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
                offset_func(exon.end if flipped else exon.start),
                offset_func(exon.start if flipped else exon.end),
                strand=new_strand,
                name=exon.name
            )
            new_exons.append(new_exon)
        new_pre_transcript = _genomic.PreTranscript(
            exons=new_exons,
            gene=new_gene,
            name=transcript.name,
            is_best_transcript=transcript.is_best_transcript
        )
        new_gene.unspliced_transcripts.append(new_pre_transcript)
        for splicing_pattern in new_pre_transcript.generate_splicing_patterns():
            new_pre_transcript.transcripts.append(
                _genomic.Transcript(new_pre_transcript, splicing_pattern)
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
