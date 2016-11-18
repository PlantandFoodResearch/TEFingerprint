#! /usr/bin/env python

import sys
import os
import argparse
from itertools import product
from multiprocessing import Pool
from tectoolkit import io
from tectoolkit.classes import ReadGroup, ReadLoci
from tectoolkit.gff import GffFeature
from tectoolkit.cluster import HUDC
from tectoolkit.cluster import FUDC


class FingerprintProgram(object):
    """Main class for the fingerprint program"""
    def __init__(self, arguments):
        if self.args.references == ['']:
            self.references = io.read_bam_references(self.args.input_bam)
        else:
            self.references = self.args.references
        """
        Init method for :class:`FingerprintProgram`.

        :param arguments: A list of commandline arguments to be parsed for the fingerprint program
        """
        self.args = self.parse_args(arguments)

    def parse_args(self, args):
        """
        Defines an argument parser to handle commandline inputs for the fingerprint program.

        :param args: A list of commandline arguments for the fingerprint program

        :return: A dictionary like object of arguments and values for the fingerprint program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('input_bam',
                            nargs=1,
                            help='A single bam file to be fingerprinted')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. ' + \
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE grouping(s) to be used. Must be exact string match(s) to start of read name')
        parser.add_argument('-s', '--strands',
                            type=str,
                            nargs='+',
                            choices=set("+-."),
                            default=['+', '-'],
                            help='Strand(s) to be analysed. Use + for forward, - for reverse and . for both')
        parser.add_argument('-e', '--eps',
                            type=int,
                            default=[100],
                            nargs='+',
                            help='Maximum allowable distance among read tips to be considered a cluster')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[5],
                            nargs=1,
                            help='Minimum allowable number of read tips to be considered a cluster')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=1,
                            help='Maximum number of cpu threads to be used')
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

    def _build_jobs(self):
        """

        :param input_bam:
        :param references:
        :param families:
        :param strands:
        :param eps:
        :param min_reads:
        :return:
        """
        return product(self.args.input_bam,
                       self.references,
                       self.args.families,
                       self.args.strands,
                       [self.args.eps],
                       self.args.min_reads)

    def _fingerprint(self, input_bam, reference, family, strand, eps, min_reads):
        """

        :param input_bam:
        :param reference:
        :param family:
        :param strand:
        :param eps:
        :param min_reads:
        :return:
        """
        fingerprint = Fingerprint(input_bam, reference, family, strand, eps, min_reads)
        for feature in fingerprint.loci_to_gff():
            print(feature)

    def run(self):
        """

        :return:
        """
        jobs = self._build_jobs()
        if self.args.threads == 1:
            for job in jobs:
                self._fingerprint(*job)
        else:
            with Pool(self.args.threads) as pool:
                pool.starmap(self._fingerprint, jobs)


class Fingerprint(object):

    def __init__(self, bam, reference, family, strand, eps, min_reads):
        self.reference = reference
        self.family = family
        self.strand = strand
        self.eps = eps
        self.min_reads = min_reads
        self.sample_name = ''
        self.source = os.path.basename(bam)
        self.reads = self._reads_from_bam(bam)
        self.loci = self._fit()

    def _reads_from_bam(self, bam):
        sam = io.read_bam_strings(bam, reference=self.reference, family=self.family, strand=self.strand)
        return ReadGroup.from_sam_strings(sam, strand=self.strand)

    def _fit(self):
        if len(self.eps) == 2:
            # use hierarchical clustering method
            max_eps, min_eps = max(self.eps), min(self.eps)
            hudc = HUDC(self.min_reads, max_eps, min_eps)
            hudc.fit(self.reads['tip'])
            return ReadLoci(hudc.loci)
        elif len(self.eps) == 1:
            # use flat clustering method
            eps = max(self.eps)
            fudc = FUDC(self.min_reads, eps)
            fudc.fit(self.reads['tip'])
            return ReadLoci(fudc.loci)
        else:
            pass

    def loci_to_gff(self):
        for start, end in self.loci:
            yield GffFeature(seqid=self.reference,
                             start=start,
                             end=end,
                             strand=self.strand,
                             ID="{0}_{1}_{2}_{3}".format(self.family, self.reference, self.strand, start),
                             Name=self.family,
                             sample=self.source)

if __name__ == '__main__':
    pass
