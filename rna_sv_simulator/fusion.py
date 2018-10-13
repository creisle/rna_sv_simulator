'''
create the fusion chromosome and re-mapping of the genes on the mutant chrs
'''
from mavis.constants import SVTYPE
from mavis.annotate import fusion as _fusion


def _mutate_translocation(reference_genome, annotations, breakpoint_pair):
    pass


def shift_gene(gene, offset):
    pass


def _mutate_continuous(reference_genome, annotations, breakpoint_pair):
    # create the original fusion transcirpt
    ft = _fusion.FusionTranscript.build(...)

    mutant_genes = []
    mutant_seq = ''

    for gene in annotations[breakpoint_pair.break1.chr]:
        if gene.end < breakpoint_pair.break1.start:
            # all genes to the left remain the same
            mutant_genes.append(gene)
        elif gene.start > breakpoint_pair.break2.end:
            # genes to the right are affected by the event
            mutant_gene = shift_gene(gene, ...) #TODO
        else:
            # genes in the middle depend on the event type


def mutate(reference_genome, annotations, breakpoint_pair):
    '''
    Returns the mutated reference genome sequences as well as the
    shifted annotations and the fusion transcript
    '''
    pass


