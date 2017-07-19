#! /usr/bin/env python

import argparse
import pysam
import subprocess
import os
import shutil
import sys
import re
from tempfile import mkdtemp


class PreProcessProgram(object):
    """"""
    def __init__(self, input_bam,
                 reference_fasta,
                 output_bam,
                 include_tails=True,
                 tail_minimum_length=38,
                 mate_element_tag='ME',
                 temp_dir=False,
                 threads=1):
        self.input_bam = input_bam
        self.reference_fasta = reference_fasta
        self.output_bam = output_bam
        self.include_tails = include_tails
        self.tail_minimum_length = tail_minimum_length
        self.mate_element_tag = mate_element_tag
        self.TEMP_DIR = temp_dir
        self.threads = str(threads)

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

    def _run_pipeline(self, SCRATCH_DIR):
        """"""

        if self.include_tails:
            jobs = ProcessForwardDanglers, ProcessReverseDanglers, ProcessForwardTailClips, ProcessReverseTailClips
        else:
            jobs = ProcessForwardDanglers, ProcessReverseDanglers

        bams = []  # used to accumulate names of bam from each job

        for job in jobs:
            # run each job and append resulting bam file to the list
            bams.append(job(self.input_bam,
                            self.reference_fasta,
                            SCRATCH_DIR).run())

        print('>>> Writing: {0}'.format(self.output_bam))
        merge_bams(self.output_bam, *bams)

        print('>>> Indexing: {0}'.format(self.output_bam))
        index_bam(self.output_bam)

    def run(self):
        """"""
        # check if samtools and bwa are available
        _check_programs_installed('bwa', 'samtools')

        # create temp dir for intermediate files unless one is supplied by user
        if self.TEMP_DIR:
            if not os.path.exists(self.TEMP_DIR):
                os.makedirs(self.TEMP_DIR)
            temp_dir = self.TEMP_DIR
            print('>>> Using directory at: {0}'.format(temp_dir))
        else:
            temp_dir = mkdtemp()
            print('>>> Creating temporary directory at: {0}'.format(temp_dir))

        # attempt running pipeline
        print('>>> Running Pipeline')
        try:
            self._run_pipeline(temp_dir)
        except:
            # remove temp dir unless it was supplied by user
            if self.TEMP_DIR:
                print('>>> Error Encountered... temporary files kept at: {0}'.format(temp_dir))
            else:
                print('>>> Error Encountered... Removing temporary directory at: {0}'.format(temp_dir))
                shutil.rmtree(temp_dir)
            # re-raise the error for debugging
            raise

        # remove temp dir unless it was supplied by user
        if self.TEMP_DIR:
            pass
        else:
            print('>>> Removing temporary directory at: {0}'.format(temp_dir))
            shutil.rmtree(temp_dir)


def _check_programs_installed(*args):
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


def index_bam(bam):
    """
    Create index files for bwa mem

    :param bam: bam file to index
    :type bam: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('samtools')
    return subprocess.call(['samtools index ' + bam], shell=True)


def merge_bams(output_bam, *args):
    input_bams = ' '.join(args)
    return subprocess.call(['samtools merge -f {0} {1}'.format(output_bam, input_bams)], shell=True)


def slice_bam(input_bam, flags, output_bam):
    _check_programs_installed('samtools')
    if isinstance(flags, list):
        flags = ' '.join(flags)
    command = 'samtools view -Su {0} {1} | samtools sort - -o {2}'.format(flags, input_bam, output_bam)
    return subprocess.call([command], shell=True)


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
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('bwa', 'samtools')
    command = ' '.join(['bwa', 'mem', '-t', threads, repeats_fasta, fastq_1, fastq_2,
                        '| samtools view -Su - | samtools sort - -o', output_bam])
    return subprocess.call([command], shell=True)


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


def parse_sam_string(string):
    """
    Parses the first 11 columns of a sam formatted string into a dictionary.

    :param string: sam formatted read
    :type string: str

    :return: dictionary of named read attributes
    :type: dict[str, str]
    """
    attributes = string.split("\t")
    return {'QNAME': attributes[0],
            'FLAG': attributes[1],
            'RNAME': attributes[2],
            'POS': attributes[3],
            'MAPQ': attributes[4],
            'CIGAR': attributes[5],
            'RNEXT': attributes[6],
            'PNEXT': attributes[7],
            'TLEN': attributes[8],
            'SEQ': attributes[9],
            'QUAL': attributes[10]}


def reads_as_fastq(read_dicts):
    for read in read_dicts:
        yield "@{0}\n{1}\n+\n{2}\n".format(read['QNAME'],
                                           read['SEQ'],
                                           read['QUAL'])


def get_mate_elements(read_dicts):
    return {read['QNAME']: read['RNAME'] for read in read_dicts}


class _ProcessDanglers(object):
    EXTRACT_DANGLER_FLAGS = ''
    TEMP_FILE_PREFIX = '_'

    def __init__(self,
                 input_bam,
                 reference,
                 directory,
                 remove_files=True,
                 bwa_threads=1,
                 element_tag='ME',
                 soft_clip_min_length=38):
        self.input_bam = input_bam
        self.reference = reference
        self.directory = directory
        self.remove_files = remove_files
        self.temp_dangler_bam = os.path.join(self.directory,
                                             self.TEMP_FILE_PREFIX +
                                             '0_danglers.bam')
        self.temp_fastq = os.path.join(self.directory,
                                       self.TEMP_FILE_PREFIX +
                                       '1_danglers.fastq')
        self.temp_mapped_danglers_bam = os.path.join(self.directory,
                                                     self.TEMP_FILE_PREFIX +
                                                     '2_mapped_danglers.bam')
        self.temp_tagged_mapped_danglers_bam = os.path.join(self.directory,
                                                            self.TEMP_FILE_PREFIX +
                                                            '3_tagged_mapped_danglers.bam')
        self.bwa_threads = bwa_threads
        self.element_tag = element_tag
        self.soft_clip_min_length = soft_clip_min_length
        self.element_tags = {}

    @staticmethod
    def _index_bam(bam):
        return subprocess.check_call(['samtools index ' + bam], shell=True)

    def _subset_bam(self):
        print('subsetting {0} to {1}'.format(self.input_bam, self.temp_dangler_bam))
        command = 'samtools view -Su {0} {1} | samtools sort - -o {2}'.format(self.EXTRACT_DANGLER_FLAGS,
                                                                              self.input_bam,
                                                                              self.temp_dangler_bam)
        code = subprocess.check_call([command], shell=True)
        self._index_bam(self.temp_dangler_bam)
        print('indexed {0}'.format(self.temp_dangler_bam))
        return code

    def _parse_reads(self, reads):
        """
        Converts a iterator pysam reads to a a generator of dicts.

        Any additional processing or filtering should be done here.

        :param reads: iterator of pysam reads
        :return: generator of dicts
        :rtype: generator[str, str]
        """
        for read in reads:
            yield {'name': read.qname,
                   'element': read.reference_name,
                   'sequence': read.seq,
                   'quality': read.qual}

    def _extract_danglers(self):
        with open(self.temp_fastq, 'a') as fq:
            reads = pysam.Samfile(self.temp_dangler_bam, "rb")
            reads = self._parse_reads(reads)
            for read in reads:

                # add read name and element name to dict
                self.element_tags[read['name']] = read['element']

                # write read to fastq
                fq.write("@{0}\n{1}\n+\n{2}\n".format(read['name'],
                                                      read['sequence'],
                                                      read['quality']))

    def _map_danglers_to_reference(self):
        command = 'bwa mem -t {0} {1} {2} | samtools view -Su - | samtools sort - -o {3}'
        command = command.format(self.bwa_threads,
                                 self.reference,
                                 self.temp_fastq,
                                 self.temp_mapped_danglers_bam)
        code = subprocess.call([command], shell=True)
        self._index_bam(self.temp_mapped_danglers_bam)
        return code

    def _tag_danglers(self):
        danglers_untagged = pysam.Samfile(self.temp_mapped_danglers_bam, "rb")
        danglers_tagged = pysam.Samfile(self.temp_tagged_mapped_danglers_bam, "wb", template=danglers_untagged)
        for read in danglers_untagged.fetch():
            read.tags += [(self.element_tag, self.element_tags[read.qname])]
            danglers_tagged.write(read)
        danglers_untagged.close()
        danglers_tagged.close()
        self._index_bam(self.temp_tagged_mapped_danglers_bam)

    def run(self):
        self._subset_bam()
        self._extract_danglers()
        self._map_danglers_to_reference()
        self._tag_danglers()
        return self.temp_tagged_mapped_danglers_bam


class ProcessForwardDanglers(_ProcessDanglers):
    EXTRACT_DANGLER_FLAGS = '-F 2304 -F 16 -F 8 -f 4'
    TEMP_FILE_PREFIX = 'forward_'

    def _parse_reads(self, reads):
        """
        Converts a iterator pysam reads to a a generator of dicts.

        Any additional processing or filtering should be done here.

        :param reads: iterator of pysam reads
        :return: generator of dicts
        :rtype: generator[str, str]
        """
        for read in reads:
            yield {'name': read.qname,
                   'element': read.reference_name,
                   'sequence': read.seq,
                   'quality': read.qual}


class ProcessReverseDanglers(_ProcessDanglers):
    EXTRACT_DANGLER_FLAGS = '-F 2304 -f 16 -F 8 -f 4'
    TEMP_FILE_PREFIX = 'reverse_'

    def _parse_reads(self, reads):
        """
        Converts a iterator pysam reads to a a generator of dicts.

        Any additional processing or filtering should be done here.

        :param reads: iterator of pysam reads
        :return: generator of dicts
        :rtype: generator[str, str]
        """
        for read in reads:
            yield {'name': read.qname,
                   'element': read.reference_name,
                   'sequence': reverse_complement(read.seq),
                   'quality': read.qual[::-1]}


class ProcessForwardTailClips(_ProcessDanglers):
    EXTRACT_DANGLER_FLAGS = '-F 2304 -f 2 -F 16 -F 4 -F 8'
    TEMP_FILE_PREFIX = 'forward_tail_clips_'

    def _parse_reads(self, reads):
        """
        Converts a iterator pysam reads to a a generator of dicts.

        Any additional processing or filtering should be done here.

        :param reads: iterator of pysam reads
        :return: generator of dicts
        :rtype: generator[str, str]
        """
        for read in reads:
            clip = re.findall(r"^[0-9]+[S]", read.cigarstring)
            # if there is a soft clipped section
            if clip:
                # get clip length
                clip = int(clip[0].strip('S'))
                if clip >= self.soft_clip_min_length:
                    # return the clipped section as a read
                    yield {'name': read.qname,
                           'element': read.reference_name,
                           'sequence': read.seq[0: clip + 1],
                           'quality': read.qual[0: clip + 1]}
                else:
                    pass
            else:
                pass


class ProcessReverseTailClips(_ProcessDanglers):
    EXTRACT_DANGLER_FLAGS = '-F 2304 -f 2 -f 16 -F 4 -F 8'
    TEMP_FILE_PREFIX = 'reverse_tail_clips_'

    def _parse_reads(self, reads):
        """
        Converts a iterator pysam reads to a a generator of dicts.

        Any additional processing or filtering should be done here.

        :param reads: iterator of pysam reads
        :return: generator of dicts
        :rtype: generator[str, str]
        """
        for read in reads:
            clip = re.findall(r"[0-9]+[S]$", read.cigarstring)
            # if there is a soft clipped section
            if clip:
                # get clip length
                clip = int(clip[0].strip('S'))
                if clip >= self.soft_clip_min_length:
                    # return the clipped section as a read
                    yield {'name': read.qname,
                           'element': read.reference_name,
                           'sequence': reverse_complement(read.seq[0: clip + 1]),
                           'quality': read.qual[0: clip + 1][::-1]}
                else:
                    pass
            else:
                pass


def extract_forward_danglers(bam):
    """
    Extracts unmapped reads with a mapped pair where the unmapped read does not have the reverse bit flag.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: list[str]
    """
    for line in pysam.view('-F', '2304', '-F', '16', '-F', '8', '-f', '4', bam).splitlines():
        read = parse_sam_string(line)
        yield read


def extract_reverse_danglers(bam):
    """
    Extracts unmapped reads with a mapped pair where the unmapped read does have the reverse bit flag.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: list[str]
    """
    for line in pysam.view('-F', '2304', '-f', '16', '-F', '8', '-f', '4', bam).splitlines():
        read = parse_sam_string(line)
        read['SEQ'] = reverse_complement(read['SEQ'])
        read['QUAL'] = read['QUAL'][::-1]
        yield read


def extract_forward_soft_clips(bam, min_length=38):
    """
    Extracts forward reads of proper mapped pairs.

    When mapping reads to repeat elements, sometimes both reads will map as a proper pair.
    If one of these reads has a soft-clipped tail, it can be used as an additional 'dangler' read.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: list[str]
    """
    for line in pysam.view('-F', '2304', '-f', '2', '-F', '16', '-F', '4', '-F', '8', bam).splitlines():
        read = parse_sam_string(line)
        # extract soft clipped tail
        clip = re.findall(r"^[0-9]+[S]", read['CIGAR'])
        if clip:
            # get clip length
            clip = int(clip[0].strip('S'))
            if clip >= min_length:
                # use soft clipped section as read sequence and quality
                read['SEQ'] = read['SEQ'][0: clip + 1]
                read['QUAL'] = read['QUAL'][0: clip + 1]
                yield read
            else:
                pass
        else:
            pass


def extract_reverse_soft_clips(bam, min_length=38):
    """
    Extracts reverse reads of proper mapped pairs.

    When mapping reads to repeat elements, sometimes both reads will map as a proper pair.
    If one of these reads has a soft-clipped tail, it can be used as an additional 'dangler' read.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: list[str]
    """
    for line in pysam.view('-F', '2304', '-f', '2', '-f', '16', '-F', '4', '-F', '8', bam).splitlines():
        read = parse_sam_string(line)
        # extract soft clipped tail
        clip = re.findall(r"[0-9]+[S]$", read['CIGAR'])
        if clip:
            # get clip length
            clip = int(clip[0].strip('S'))
            if clip >= min_length:
                # use soft clipped section as read sequence and quality
                read['SEQ'] = reverse_complement(read['SEQ'][-clip:])
                read['QUAL'] = read['QUAL'][-clip:][::-1]
                yield read
            else:
                pass
        else:
            pass


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
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('bwa', 'samtools')
    command = ' '.join(['bwa', 'mem', '-t', threads, reference_fasta, fastq,
                        '| samtools view -Su - | samtools sort - -o', output_bam])
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
