#! /usr/bin/env python

import numpy as np
from tectoolkit import bam_io


class ReadGroup(object):
    """
    A collection of mapped SAM read positions. This class does not contain read sequences.

    Each read is represented as a single element in a numpy array with the following fields:
     - tip: np.int64
     - tail: np.int64
     - strand: np.str_, 1
     - name: np.str_, 256
    Where 'tip' and 'tail' are the coordinates (1 based indexing) of the read mapped to its reference,
    'strand' is a single character ('-' or '+') indicating the strand the read mapped to
    and 'name' is a string of up to 256 characters containing the read name.
    """
    DTYPE_READ = np.dtype([('tip', np.int64),
                           ('tail', np.int64),
                           ('strand', np.str_, 1),
                           ('name', np.str_, 254)])

    def __init__(self, reads=None):
        """
        Init method for :class:`ReadGroup`.

        :param reads: A numpy array.
        :type reads: :class:`numpy.ndarray`[(int, int, str, str)]
        """
        if reads is None:
            reads = []
        self.reads = np.array(reads, dtype=ReadGroup.DTYPE_READ, copy=True)

    def __iter__(self):
        """
        Iter method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :return: An iterable of mapped SAM read positions and names
        :rtype: generator[(int, int, str, str)]
        """
        for read in self.reads:
            yield read

    def __getitem__(self, item):
        """
        Getitem method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :param item:
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`ReadGroup`.DTYPE_READ
        :rtype: :class:`numpy.ndarray`[(int, int, str, str)]
        """
        return self.reads[item]

    def __len__(self):
        """
        Len method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :return: Number of reads in group
        :rtype: int
        """
        return len(self.reads)

    def sort(self, order='tip'):
        """
        Sort reads in place by field(s).

        :param order: A valid field or list of fields in :class:`ReadGroup`, defaults to 'tip'
        :type order: str | list[str]
        """
        self.reads.sort(order=order)

    def strand(self):
        """
        Identifies the consensus strand of all reads in the group.

        :return: '+', '-' or '.' respectively if all reads are on forwards, reverse or combination of strands
        :rtype: str
        """
        variation = set(self.reads['strand'])
        if len(variation) > 1:
            return '.'
        else:
            return variation.pop()

    def subset_by_locus(self, start, stop, margin=0, end='tip'):
        """
        Returns a new ReadGroup object containing (the specified end of) all reads within specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: The read end that must fall within the bounds, must be 'tip' or 'tail', defaults to 'tip'
        :type end: str

        :return: The subset of reads that fall within the specified bounds
        :rtype: :class:`ReadGroup`
        """
        assert end in {'tip', 'tail', 'both'}
        start -= margin
        stop += margin
        if end == 'both':
            # 'tip' may be smaller or larger than 'tail'
            reads = self.reads
            reads = reads[np.logical_and(reads['tip'] >= start, reads['tip'] <= stop)]
            reads = reads[np.logical_and(reads['tail'] >= start, reads['tail'] <= stop)]
        else:
            reads = self.reads[np.logical_and(self.reads[end] >= start, self.reads[end] <= stop)]
        return ReadGroup(reads)

    @classmethod
    def _parse_sam_strings(cls, strings, strand=None):
        """
        Parses a collection of SAM formatted strings into a tuple generator.

        :param strings: A collection of SAM formatted strings
        :type strings: iterable[str]
        :param strand: Strand ('+' or '-') of all reads (if known)
        :type strand: str

        :return: An iterable of mapped SAM read positions and names
        :rtype: generator[(int, int, str, str)]
        """
        def _parse_sam_string(string, strand):
            attr = string.split("\t")
            name = str(attr[0])
            start = int(attr[3])
            length = len(attr[9])
            end = start + length - 1  # 1 based indexing used in SAM format
            if strand is None:
                strand = bam_io.flag_orientation(int(attr[1]))
            if strand == '+':
                tip = end
                tail = start
                return tip, tail, strand, name
            elif strand == '-':
                tip = start
                tail = end
                return tip, tail, strand, name

        assert strand in ['+', '-', None]
        reads = (_parse_sam_string(string, strand) for string in strings)
        return reads

    @classmethod
    def _from_sam_strings(cls, strings, strand=None):
        """
        Construct an instance of :class:`ReadGroup` from an iterable of SAM formatted strings.

        :param strings: An iterable of SAM formatted strings
        :type strings: iter[str]
        :param strand: Strand ('+' or '-') of all reads (if known)
        :type strand: str

        :return: An instance of :class:`ReadGroup`
        :rtype: :class:`ReadGroup`
        """
        generator = cls._parse_sam_strings(strings, strand=strand)
        return cls.from_iter(generator)

    @classmethod
    def from_iter(cls, iterable):
        """
        Create an instance of :class:`ReadGroup` from an iterable.

        :param iterable: an iterable of read positions
        :type iterable: iter[(int, int, str, str)]

        :return: An instance of :class:`ReadGroup`
        :rtype: :class:`ReadGroup`
        """
        reads = np.fromiter(iterable, dtype=ReadGroup.DTYPE_READ)
        reads.sort(order=('tip', 'tail'))
        return ReadGroup(reads)

    @classmethod
    def from_bam(cls, bam, reference, family, strand):
        """
        Reads subset of reads from bam file that are of the targeted reference, family and strand.

        :param bam: A bam file path
        :type bam: str
        :param reference: Target reference of reads
        :type reference: str
        :param family: Target family of reads
        :type family: str
        :param strand: Target strand of reads ('+' or '-')
        :type strand: str

        :return: Subset of bam reads
        :rtype: :class:`ReadGroup`
        """
        assert strand in ('+', '-')
        sam = bam_io.read_bam_strings(bam, reference=reference, family=family, strand=strand)
        return ReadGroup._from_sam_strings(sam, strand=strand)

if __name__ == '__main__':
    pass