#! /usr/bin/env python

import pysam
import numpy as np

locus = np.dtype([('start', np.int64),
                  ('stop', np.int64)])

sam_read = np.dtype([('tip', np.int64),
                     ('tail', np.int64),
                     ('length', np.int64),
                     ('reverse', np.bool)])


class GffFeature(object):
    """"""
    def __init__(self, seqid, source='.', ftype='.', start='.', end='.', score='.', strand='.', phase='.', **kwargs):
        self.seqid = seqid
        self.source = source
        self.type = ftype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = kwargs

    def attribute_names(self):
        return set(self.attributes.keys())

    def __parse_attributes(self, attributes):
        if attributes is not None:
            return ';'.join(tuple('{0}={1}'.format(key, value) for key, value in attributes.items()))
        else:
            return'.'

    def __str__(self):
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(self.seqid,
                               self.source,
                               self.type,
                               self.start,
                               self.end,
                               self.score,
                               self.strand,
                               self.phase,
                               self.__parse_attributes(self.attributes))


def strand2flag(strand):
    """

    :param strand:
    :return:
    """
    if strand == '+':
        return ('-F', '20')
    elif strand == '-':
        return ('-f', '16')
    elif strand == '.':
        return ('-F', '4')
    else:
        pass  # throw error


def read_sam_strings(input_bam, reference, family, strand):
    """

    :param reference:
    :param family:
    :param strand:
    :return:
    """
    flag = strand2flag(strand)
    sam_strings = (string for string in pysam.view(flag[0], flag[1], input_bam, reference, family).splitlines())
    return sam_strings


def read_bam_references(input_bam):
    bam = pysam.AlignmentFile(input_bam, 'rb')
    references = bam.references
    bam.close()
    return references


def parse_sam_strings(sam_strings, strand):
    """

    :param sam_strings:
    :param strand:
    :return:
    """
    def parse_sam_string(sam_string, strand):
        attr = sam_string.split("\t")
        start = int(attr[3])
        length = len(attr[9])
        end = start+length
        if strand == '+':
            reverse = True
            tip = end
            tail = start
            return tip, tail, length, reverse
        elif strand == '-':
            reverse = False
            tip = start
            tail = end
            return tip, tail, length, reverse

    reads = (parse_sam_string(string, strand) for string in sam_strings)
    reads = np.fromiter(reads, dtype=sam_read)
    reads.sort(order=('tip', 'tail'))
    return reads


def read_sam_reads(input_bam, reference, family, strand):
    """

    :param reference:
    :param family:
    :param strand:
    :return:
    """
    sam_strings = read_sam_strings(input_bam, reference, family, strand)
    sam_reads = parse_sam_strings(sam_strings, strand)
    return sam_reads


def simple_subcluster(reads, minpts, eps):
    """

    :param reads:
    :param minpts:
    :param eps:
    :return:
    """
    reads.sort(order=('tip', 'tail'))
    offset = minpts - 1
    upper = reads['tip'][offset:]
    lower = reads['tip'][:-offset]
    diff = upper - lower
    dense = diff <= eps
    lower = lower[dense]
    upper = upper[dense]
    loci = ((lower[i], upper[i]) for i in range(len(lower)))
    loci = np.fromiter(loci, dtype=locus)
    return loci


def merge_clusters(clusters):
    """

    :param clusters:
    :return:
    """
    def merge(sub_clusters):
        sub_clusters.sort(order=('start', 'stop'))
        start = sub_clusters['start'][0]
        stop = sub_clusters['stop'][0]
        for i in range(1, len(sub_clusters)):
            if sub_clusters['start'][i] <= stop:
                stop = sub_clusters['stop'][i]
            else:
                yield start, stop
                start = sub_clusters['start'][i]
                stop = sub_clusters['stop'][i]
        yield start, stop
    return np.fromiter(merge(clusters), dtype=locus)


def simple_cluster(reads, minpts, eps):
    """

    :param points:
    :param minpts:
    :param eps:
    :return:
    """
    subclusters = simple_subcluster(reads, minpts, eps)
    return merge_clusters(subclusters)


def resample_reads(reads, n=None, reps=1):
    """

    :param reads:
    :param n:
    :param reps:
    :return:
    """
    length = len(reads)
    if n is None:
        n = length

    def resample(sample, n):
        indexes = np.random.choice(len(sample), n)
        new_sample = sample[indexes]
        new_sample.sort(order='tip')
        return new_sample
    return (resample(reads, n) for _ in range(reps))


def bootstrapped_simple_cluster(reads, minpts, eps, reps, p=95):
    """

    :param reads:
    :param minpts:
    :param eps:
    :param reps:
    :param support:
    :return:
    """
    sample_sets = resample_reads(reads, reps=reps)
    replicates = (simple_cluster(sample, minpts, eps) for sample in sample_sets)
    depth = np.zeros(max(reads['tip']))
    for clusters in replicates:
        for start, stop in clusters:
            depth[start:(stop+1)] += 1
    depth = (depth/reps) * 100
    supported_points = np.where(depth >= p)[0]
    supported_points.sort()
    supported_clusters = simple_cluster(supported_points, minpts, eps)
    return supported_clusters


def reads_in_locus(reads, locus, margin=0):
    """

    :param cluster:
    :param points:
    :return:
    """
    start, stop = locus
    subset = reads[(reads['tip'] >= start - margin) & (reads['tip'] <= stop + margin)]
    return subset


def cluster_depth(reads, cluster):
    """

    :param reads:
    :param cluster:
    :return:
    """
    start, stop = cluster
    depth = np.zeros(stop - (start-1), dtype=np.int)
    tips = reads['tip'] - (start)
    for tip in tips:
        depth[tip] += 1
    return depth


if __name__ == '__main__':
    pass