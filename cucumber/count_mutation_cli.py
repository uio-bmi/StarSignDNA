import bionumpy as bnp
import numpy as np
from bionumpy.variants import count_mutation_types_genomic, count_mutation_types, MutationTypeEncoding
from bionumpy.io.delimited_buffers import PhasedVCFMatrixBuffer
from bionumpy.io.matrix_dump import matrix_to_csv, read_matrix
from bionumpy import Genome
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename: str, fasta_filename: str, out_filename: str = None, flank: int = 1, has_numeric_chromosomes=True, genotyped=False):
    buffer_type = PhasedVCFMatrixBuffer if genotyped else None
    genome = Genome.from_file(fasta_filename)
    variants = np.concatenate(list(bnp.open(vcf_filename, buffer_type=buffer_type).read_chunks()))
    variants = genome.get_locations(variants, has_numeric_chromosomes=has_numeric_chromosomes)
    counts = count_mutation_types_genomic(variants, genome.read_sequence(), genotyped=genotyped)
    if out_filename is not None:
        output = matrix_to_csv(counts.counts, header=counts.alphabet)
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(counts)
    return counts

def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa", has_numeric_chromosomes=False)


def pipeline(vcf_filename: str, fasta_filename: str, signature_filename, has_numeric_chromosomes=True, genotyped=False):
    counts = main(vcf_filename, fasta_filename, has_numeric_chromosomes=has_numeric_chromosomes, genotyped=genotyped)
    signature_matrix = read_matrix(signature_filename)
    count_matrix = signature_matrix.__class__(np.atleast_2d(counts.counts), col_names=counts.alphabet)
    print(nmf(signature_matrix, count_matrix))


if __name__ == "__main__":
    import typer
    typer.run(pipeline)
    # run_as_commandline(main)

