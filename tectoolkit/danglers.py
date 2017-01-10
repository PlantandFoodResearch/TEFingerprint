#! /usr/bin/env python

import pysam
import subprocess
import os
from tempfile import mkdtemp


class PreProcessProgram(object):
    """"""

    def __init__(self, fastq_1, fastq_2, reference_fasta, repeats_fasta, output_bam, threads=1):
        self.fastq_1 = fastq_1
        self.fastq_2 = fastq_2
        self.reference_fasta = reference_fasta
        self.repeats_fasta = repeats_fasta
        self.output_bam = output_bam
        self.threads = threads

    def run(self):
        """"""

        # create temp dir for intermediate files
        temp_dir = mkdtemp()

        # map reads to repeats and store as temp file
        temp_bam_1 = os.path.join(temp_dir, 'temp_bam_1')
        map_pairs_to_repeat_elements(self.fastq_1,
                                     self.fastq_2,
                                     self.repeats_fasta,
                                     temp_bam_1,
                                     self.threads)

        # extract danglers from bam
        dangler_strings = extract_danglers(temp_bam_1)

        # convert dangler strings to temp fastq
        temp_fastq_danglers = os.path.join(temp_dir, 'temp_fastq_danglers')
        sam_strings_to_fastq(dangler_strings,
                             temp_fastq_danglers)

        # map danglers to reference
        temp_bam_2 = os.path.join(temp_dir, 'temp_bam_2')
        map_danglers_to_reference(temp_fastq_danglers,
                                  self.reference_fasta,
                                  temp_bam_2,
                                  self.threads)

        # tag danglers and write output file
        tag_danglers(temp_bam_2, dangler_strings, self.output_bam)

        # remove temp dir
        os.rmdir(temp_dir)


def map_pairs_to_repeat_elements(fastq_1, fastq_2, repeats_fasta, output_bam, threads=1):
    """
    Maps paired end reads to repeat-element reference and writes bam file.

    :param fastq_1: paired end reads file 1
    :type fastq_1: str
    :param fastq_2: paired end reads file 2
    :type fastq_2: str
    :param repeats_fasta: repeat-element reference
    :type repeats_fasta: str
    :param output_bam: bam file of paired end reads aligned to repeat-element reference
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: int
    """
    subprocess.call(['bwa', 'mem', '-t', threads, repeats_fasta, fastq_1, fastq_2,
                     '| samtools view -Su - | samtools sort - -o', output_bam],
                    shell=True)


def extract_danglers(bam):
    """
    Extracts unmapped reads with a mapped pair.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: str
    """
    return pysam.samtools.view('-F', '0x800', '-F', '0x100', '-F', '8', '-f', '4', bam)


def sam_strings_to_fastq(sam_strings, output_fastq):
    """
    Extracts unmapped reads with a mapped pair and writes them to a fastq file.

    :param sam_strings: sam formatted strings
    :type sam_strings: str
    :param output_fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type output_fastq: str
    """
    sam_attributes = iter(line.split('\t') for line in sam_strings.splitlines())
    fastq_lines = ("@{0}\n{1}\n+\n{2}\n".format(items[0], items[9], items[10]) for items in sam_attributes)
    with open(output_fastq, 'w') as f:
        for line in fastq_lines:
            f.write(line)


def map_danglers_to_reference(fastq, reference_fasta, output_bam, threads):
    """
    Maps dangler reads (non-mapped reads with pair mapping to repeat-element) to reference genome and writes a bam file.

    :param fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type fastq: str
    :param reference_fasta: fasta file containing reference genome
    :type reference_fasta: str
    :param output_bam: bam file of dangler reads aligned to reference genome
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: int
    """
    subprocess.call(['bwa', 'mem', '-t', threads, reference_fasta, fastq,
                     '| samtools view -Su - | samtools sort - -o', output_bam],
                    shell=True)


def tag_danglers(dangler_bam, dangler_strings, output_bam):
    """
    Tags mapped danglers with the id of the element which their mate mapped to and writes to new bam.

    :param dangler_bam: bam file of dangler reads aligned to reference genome
    :type dangler_bam: str
    :param dangler_strings: sam formatted strings of non-mapped reads with pair mapping to repeat-element
    :type dangler_strings: str
    :param output_bam: bam file of dangler reads aligned to reference genome and tagged with mate element
    :type output_bam: str
    """
    sam_attributes = iter(line.split('\t') for line in dangler_strings.splitlines())
    read_mate_elements = {attrs[0]: attrs[2] for attrs in sam_attributes}
    danglers_untagged = pysam.Samfile(dangler_bam, "rb")
    danglers_tagged = pysam.Samfile(output_bam, "wb", template=danglers_untagged)
    for read in danglers_untagged.fetch():
        read.tags += [('ME', read_mate_elements[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    pass
