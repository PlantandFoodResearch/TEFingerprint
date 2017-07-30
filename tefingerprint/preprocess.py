#! /usr/bin/env python

import argparse
import pysam
import subprocess
import os
import shutil
import sys
from tempfile import mkdtemp


class PreProcessProgram(object):
    """"""
    def __init__(self, input_bam,
                 reference_fasta,
                 output_bam,
                 include_tails=True,
                 tail_minimum_length=38,
                 mate_element_tag='ME',
                 temp_dir=None,
                 threads=1):
        self.input_bam = input_bam
        self.reference_fasta = reference_fasta
        self.output_bam = output_bam
        self.include_tails = include_tails
        self.tail_minimum_length = tail_minimum_length
        self.mate_element_tag = mate_element_tag
        self.temp_dir = temp_dir
        self.threads = str(threads)
        if self.temp_dir:
            self.keep_temp_dir = True
        else:
            self.keep_temp_dir = False

    @staticmethod
    def _cli(args):
        """
        Argument parser to handle commandline inputs for the preprocessing program.

        :param args: List of commandline arguments for the preprocess program
        :type args: []

        :return: Dictionary like object of arguments and values
        """
        parser = argparse.ArgumentParser('Extract reads whose pair was mapped to a repeat element, '
                                         'map them to a reference genome and tag them with their mates element name.')
        parser.add_argument('bam',
                            type=str,
                            nargs=1,
                            help='A bam file containing paired end reads mapped to a library of repeat elements.')
        parser.add_argument('-r', '--reference',
                            type=str,
                            nargs=1,
                            help="Reference genome in fasta format indexed with bwa.")
        parser.add_argument('-o', '--output',
                            type=str,
                            nargs=1,
                            help="Name of output bam file.")
        parser.set_defaults(include_tails=True)
        parser.add_argument('--include-tails',
                            dest='include_tails',
                            action='store_true',
                            help="Include soft-clipped tails of pairs properly mapping to a repeat element (default).")
        parser.add_argument('--exclude-tails',
                            dest='include_tails',
                            action='store_false',
                            help="Don't include soft-clipped tails.")
        parser.add_argument('--tail-minimum-length',
                            nargs=1,
                            type=int,
                            default=[38],
                            help="Minimum required length for inclusion of soft-clipped tails.")
        parser.add_argument('--mate-element-tag',
                            type=str,
                            nargs=1,
                            default=['ME'],
                            help='Tag used in bam file to indicate the element type mate read.')
        parser.add_argument('--tempdir',
                            type=str,
                            nargs=1,
                            default=[False],
                            help="Optional directory to store temp files in (files will not be removed).")
        parser.add_argument('-t', '--threads',
                            type=int,
                            nargs=1,
                            default=[1],
                            help='Maximum number of cpu threads to be used.')
        return parser.parse_args(args)

    @staticmethod
    def from_cli(args):
        arguments = PreProcessProgram._cli(args)
        return PreProcessProgram(arguments.bam[0],
                                 arguments.reference[0],
                                 arguments.output[0],
                                 include_tails=arguments.include_tails,
                                 tail_minimum_length=arguments.tail_minimum_length[0],
                                 mate_element_tag=arguments.mate_element_tag[0],
                                 temp_dir=arguments.tempdir[0],
                                 threads=arguments.threads[0])

    def _run_pipeline(self):
        """"""
        temp_fastq = self.temp_dir + '/0.fastq'
        temp_bam = self.temp_dir + '/1.bam'

        print('>>> Extracting reads from : {0}'.format(self.input_bam))
        informative_reads = extract_informative_reads(self.input_bam,
                                                      include_tails=self.include_tails,
                                                      minimum_tail=self.tail_minimum_length)

        print('>>> Writing reads to : {0}'.format(temp_fastq))
        element_tags = capture_elements_and_write_fastq(informative_reads, temp_fastq)

        print('>>> Mapping reads to {0} in {1}'.format(self.reference_fasta, temp_bam))
        map_danglers_to_reference(temp_fastq, self.reference_fasta, temp_bam, self.threads)

        print('>>> Indexing: {0}'.format(temp_bam))
        pysam.index(temp_bam)

        print('>>> Tagging reads in: {0}'.format(self.output_bam))
        tag_danglers(temp_bam, element_tags, self.output_bam, self.mate_element_tag)

        print('>>> Indexing: {0}'.format(self.output_bam))
        pysam.index(self.output_bam)

    def run(self):
        """"""
        # check if samtools and bwa are available
        check_programs_installed('bwa', 'samtools')

        # create temp dir for intermediate files unless one is supplied by user
        if self.keep_temp_dir:
            print('>>> Using directory at: {0}'.format(self.temp_dir))
            if not os.path.exists(self.temp_dir):
                os.makedirs(self.temp_dir)
        else:
            print('>>> Creating temporary directory at: {0}'.format(self.temp_dir))
            self.temp_dir = mkdtemp()

        # attempt running pipeline
        print('>>> Running Pipeline')
        try:
            self._run_pipeline()
        except:
            # remove temp dir unless it was supplied by user
            if self.temp_dir:
                print('>>> Error Encountered... temporary files kept at: {0}'.format(self.temp_dir))
            else:
                print('>>> Error Encountered... Removing temporary directory at: {0}'.format(self.temp_dir))
                shutil.rmtree(self.temp_dir)
            # re-raise the error for debugging
            raise

        # remove temp dir unless it was supplied by user
        if self.keep_temp_dir:
            pass
        else:
            print('>>> Removing temporary directory at: {0}'.format(self.temp_dir))
            shutil.rmtree(self.temp_dir)


def extract_informative_reads(bam, include_tails=True, minimum_tail=38):
    """
    Extracts informative reads from a bam file of paired end reads aligned to a library of transposons.

    :param bam: path to a bam file of paired end reads aligned to a library of transposons
    :type bam: str
    :param include_tails: whether or not to include soft clipped tails as an additional source of information
    :type include_tails: bool
    :param minimum_tail: minimum length of tail soft clips to include
    :type minimum_tail: int

    :return: parsed and filtered reads relating information about location of transposon insertions
    :rtype: generator[dict[str: str]]
    """
    alignment = pysam.AlignmentFile(bam, 'r')
    for read in alignment:
        if any((read.is_secondary, read.is_supplementary, read.mate_is_unmapped)):
            # read is bad - for our cases mate should always be mapped
            pass
        else:
            if read.is_unmapped:
                # read is a dangler

                if read.is_reverse:
                    # read is a reverse dangler
                    yield {'name': read.qname,
                           'element': read.next_reference_name,  # will be identical to read.reference_name
                           'sequence': reverse_complement(read.seq),
                           'quality': read.qual[::-1]}

                else:
                    # read is a forward dangler
                    yield {'name': read.qname,
                           'element': read.next_reference_name,  # will be identical to read.reference_name
                           'sequence': read.seq,
                           'quality': read.qual}

            elif include_tails and read.is_proper_pair:
                # potential tail clip

                if read.is_reverse:
                    # potential reverse clip
                    tail = read.cigar[-1]
                    if tail[0] == 4 and tail[1] >= minimum_tail:
                        # tail is soft clipped and of sufficient length
                        clip = tail[1]
                        # return the clipped section as a read
                        yield {'name': read.qname,
                               'element': read.reference_name,
                               'sequence': reverse_complement(read.seq[-clip:]),
                               'quality': read.qual[-clip:][::-1]}
                    else:
                        # not a soft clipped section or not  sufficient length
                        pass

                else:
                    # potential forward clip
                    tail = read.cigar[0]
                    if tail[0] == 4 and tail[1] >= minimum_tail:
                        # tail is soft clipped and of sufficient length
                        clip = tail[1]
                        # return the clipped section as a read
                        yield {'name': read.qname,
                               'element': read.reference_name,
                               'sequence': read.seq[0: clip + 1],
                               'quality': read.qual[0: clip + 1]}
                    else:
                        # not a soft clipped section or not  sufficient length
                        pass
            else:
                # neither dangler nor tail clip
                pass


def capture_elements_and_write_fastq(reads, fastq):
    """
    Captures writes parsed reads to a fastq file and returns a dictionary of read-element pairs.

    :param reads: sequence of dictionaries of reads with fields 'name', 'element', sequence' and 'quality'
    :type reads: iterator[dict[str: str]]
    :param fastq: path to fastq to write
    :type fastq: str

    :return: dictionary of read-element pairs
    :rtype: dict[str, str]
    """
    element_tags = {}

    with open(fastq, 'w') as fq:
        for read in reads:

            # add read name and element name to dict
            element_tags[read['name']] = read['element']

            # write read to fastq
            fq.write("@{0}\n{1}\n+\n{2}\n".format(read['name'],
                                                  read['sequence'],
                                                  read['quality']))
    return element_tags


def reverse_complement(sequence):
    """
    Computes the reverse compliment of a nucleotide string.

    :param sequence: capitalised string of nucleotides {'A', 'T', 'G', 'C', 'N'}
    :type sequence: str

    :return: capitalised reverse compliment string of nucleotides
    :rtype: str
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(sequence.upper()))


def check_programs_installed(*args):
    """
    Check if a program exists on the users path.

    :param args: names of programs
    :type args: str

    :return: true if all of the program are on the users path
    :rtype: bool
    """
    for program in args:
        if shutil.which(program) is None:
            raise EnvironmentError("could not find '{0}' in $PATH: {1}".format(program, os.environ['PATH']))
    return True


def map_danglers_to_reference(fastq, reference_fasta, output_bam, threads=1):
    """
    Maps dangler reads (non-mapped reads with pair mapping to repeat-element) to reference genome and writes a bam file.

    :param fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type fastq: str
    :param reference_fasta: fasta file containing reference genome
    :type reference_fasta: str
    :param output_bam: bam file of dangler reads aligned to reference genome
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    check_programs_installed('bwa', 'samtools')
    command = 'bwa mem -t {threads} {reference} {fastq} | samtools view -Su - | samtools sort - -o {bam}'
    command = command.format(threads=threads, reference=reference_fasta, fastq=fastq, bam=output_bam)
    return subprocess.call([command], shell=True)


def tag_danglers(dangler_bam, mate_element_dict, output_bam, tag):
    """
    Tags mapped danglers with the id of the element which their mate mapped to and writes to new bam.

    :param dangler_bam: bam file of dangler reads aligned to reference genome
    :type dangler_bam: str
    :param mate_element_dict: dict of read-name, mate-element pairs
    :type mate_element_dict: dict[str, str]
    :param output_bam: bam file of dangler reads aligned to reference genome and tagged with mate element
    :type output_bam: str
    """
    danglers_untagged = pysam.Samfile(dangler_bam, "rb")
    danglers_tagged = pysam.Samfile(output_bam, "wb", template=danglers_untagged)
    for read in danglers_untagged.fetch():
        read.tags += [(tag, mate_element_dict[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    PreProcessProgram.from_cli(sys.argv)
