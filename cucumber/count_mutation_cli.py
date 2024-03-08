import bionumpy as bnp
import numpy as np
from bionumpy.variants import count_mutation_types_genomic
from bionumpy.io.vcf_buffers import PhasedVCFMatrixBuffer
from bionumpy.io.matrix_dump import matrix_to_csv
from bionumpy import Genome
import logging
logging.basicConfig(level="INFO")


def count_mutation(vcf_filename: str, fasta_filename: str, out_filename: str = None, numeric_chromosomes: bool = True, genotyped: bool= True):
    """
    Count mutation types in a VCF file.

    Parameters
    ----------
    vcf_filename: str
        VCF file of the variants (possibly many samples)
    fasta_filename: str
        Fasta file of the reference genome
    out_filename: str
    numeric_chromosomes: bool
        True if chromosome names in vcf are '1', '2', '3'. False if 'chr1', 'chr2', 'chr3'
    genotyped: bool
        True if the VCF file has genotype information for many samples

    """
    buffer_type = PhasedVCFMatrixBuffer if genotyped else None
    genome = Genome.from_file(fasta_filename)
    variants = np.concatenate(list(bnp.open(vcf_filename, buffer_type=buffer_type).read_chunks()))
    variants = genome.get_locations(variants, has_numeric_chromosomes=numeric_chromosomes)
    counts = count_mutation_types_genomic(variants, genome.read_sequence(), genotyped=genotyped)
    if out_filename is not None:
        # output = matrix_to_csv(counts.counts, header=counts.alphabet)
        output = matrix_to_csv(counts.counts.astype(int), header=counts.alphabet, sep='\t')
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(counts)
    return counts
