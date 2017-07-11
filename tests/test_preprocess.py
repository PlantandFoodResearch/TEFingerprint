#! /usr/bin/env python
import os
import pytest
import subprocess
from tefingerprint import preprocess


@pytest.mark.parametrize('query,answer',
                         [('AGATC', 'GATCT'),
                          ('aGtcn', 'NGACT')])
def test_reverse_complement(query, answer):
    assert preprocess.reverse_complement(query) == answer


def test_parse_sam_string():
    query = 'read008\t0\tchr1\t2784\t60\t58M\t*\t0\t0\tCAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTACTATTTCTTGGT\t2222222222222222222222222222222222222222222222222222222222\tNA:Z:foobar'
    answer = {'QNAME': 'read008',
              'FLAG': '0',
              'RNAME': 'chr1',
              'POS': '2784',
              'MAPQ': '60',
              'CIGAR': '58M',
              'RNEXT': '*',
              'PNEXT': '0',
              'TLEN': '0',
              'SEQ': 'CAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTACTATTTCTTGGT',
              'QUAL': '2222222222222222222222222222222222222222222222222222222222'}
    assert preprocess.parse_sam_string(query) == answer


def test_foward_danglers_as_fastq():
    query = ['read007\t0\tchr1\t2779\t60\t63M\t*\t0\t0\tTCCGGCAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTTCTTTATCTTGGT\t222222222222222222222222222222222222222222222222222222222222222',
             'read008\t0\tchr1\t2784\t60\t58M\t*\t0\t0\tCAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTACTATTTCTTGGT\t2222222222222222222222222222222222222222222222222222222222',
             'read009\t0\tchr1\t2796\t60\t46M\t*\t0\t0\tATAAGTTTTCTTTCATTTTTTTATAATTGTTTTTCTTTTTCTTGGT\t2222222222222222222222222222222222222222222222']

    answer = ["@read007\nTCCGGCAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTTCTTTATCTTGGT\n+\n222222222222222222222222222222222222222222222222222222222222222",
              "@read008\nCAAGCACCTTTCATAAGTTTTCTTTCATTTTTTTATAATTGTTTTACTATTTCTTGGT\n+\n2222222222222222222222222222222222222222222222222222222222",
              "@read009\nATAAGTTTTCTTTCATTTTTTTATAATTGTTTTTCTTTTTCTTGGT\n+\n2222222222222222222222222222222222222222222222"]

    assert preprocess.forward_danglers_as_fastq(query) == answer


def test_reverse_danglers_as_fastq():
    query = ['read015\t16\tchr1\t3217\t60\t100M\t*\t0\t0\tTCTCCTTATTTAGTAGAGACTCGTTTTTAAAGATTTAGAGGGATGTTATAGTCTTTTACCGTACTTTCCCAATAGATAACATGAACCTCAGATTCGACTT\t1112222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222',
             'read016\t16\tchr1\t3226\t60\t100M\t*\t0\t0\tTTAGTAGAGACTCGTTTTTAAAGATTTAGAGGGATGTTATAGTCTTTTACCGTACTTTCCCAATAGATAACTTGAACCTCAGATTCGACTTGGATTTTCC\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222',
             'read017\t16\tchr1\t3246\t60\t100M\t*\t0\t0\tAAGATTTAGAGGGATGTTATAGTCTTTTACCGTACTTTCCCAATAGATAACATGAACCTCAGATTCGACTTTGATTTTCCACGTATTGTTTTTTCCAAGT\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222']

    answer = ["@read015\nAAGTCGAATCTGAGGTTCATGTTATCTATTGGGAAAGTACGGTAAAAGACTATAACATCCCTCTAAATCTTTAAAAACGAGTCTCTACTAAATAAGGAGA\n+\n2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222111",
              "@read016\nGGAAAATCCAAGTCGAATCTGAGGTTCAAGTTATCTATTGGGAAAGTACGGTAAAAGACTATAACATCCCTCTAAATCTTTAAAAACGAGTCTCTACTAA\n+\n2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222",
              "@read017\nACTTGGAAAAAACAATACGTGGAAAATCAAAGTCGAATCTGAGGTTCATGTTATCTATTGGGAAAGTACGGTAAAAGACTATAACATCCCTCTAAATCTT\n+\n2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"]

    assert preprocess.reverse_danglers_as_fastq(query) == answer


def test_forward_soft_clipped_tails_as_fastq():
    query = ['read001\t163\tGypsy_Gypsy26\t8339\t0\t23S61M1I15M\t=\t8783\t544\tATGTTTTTTTTTAAATAAAGTTATTGTGGACCCCGCATTTCGGCTCATGCGTTTCCCACTCGATGGCGAGCTGAAATTTTATTTAGAACAATTGATTTTT\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tNM:i:5\tMD:Z:49C0G1T11A11\tAS:i:49\tXS:i:52',
             'read002\t99\tCopia_Copia84\t6700\t0\t43S57M\t=\t6976\t376\tTTTTATGATTCCTTGACATGAAGCAACTTACATGCCCTGTAACTGAAAGGTTATAAAATAAATTTATAACACATACAGCGGAAAATTGATCTATTTTGCA\t1222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tNM:i:0\tMD:Z:5\tAS:i:57\tXS:i:57',
             'read003\t163\tCopia_Copia84\t6700\t40\t74S26M\t=\t7030\t430\tATTTTGGAAATTAAATTTAAAAGGAATTGGATTTTATGATTCCTTGACCTGAAGCAACTTAGATGCCCTGTAACTGAAAGGTTATAAAATAAATTTATAA\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tNM:i:0\tMD:Z:26\tAS:i:26\tXS:i:26']

    answer = ['@read002\nTTTTATGATTCCTTGACATGAAGCAACTTACATGCCCTGTAACT\n+\n12222222222222222222222222222222222222222222',
              '@read003\nATTTTGGAAATTAAATTTAAAAGGAATTGGATTTTATGATTCCTTGACCTGAAGCAACTTAGATGCCCTGTAACT\n+\n222222222222222222222222222222222222222222222222222222222222222222222222222']

    assert preprocess.forward_soft_clipped_tails_as_fastq(query) == answer


def test_reverse_soft_clipped_tails_as_fastq():
    query = ['read001\t83\tMULE_MUDRAVI2\t1093\t48\t35S24M41S\t=\t617\t-500\tGTTTTAACTTTAAGTTTAAATTCATGATTGATGACACGGTTTTGATTGTTTTTAAAAAGGAAAACCGTAATTGGTAACTACGATTTTATGACATGATAAT\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222221\tNM:i:0\tMD:Z:24\tAS:i:24\tXS:i:19',
             'read002\t147\tMULE_MUDRAVI2\t1093\t12\t57S28M15S\t=\t725\t-396\tAAAGTCGTACTTGGTGACTACGCTTTTAACTTTAAGTTTAAAATAATGATTGATGACACGGTTTTGATTGTTTTTAAAAAGTAAAACCGTAATTGGTAAC\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tNM:i:0\tMD:Z:28\tAS:i:28\tXS:i:29',
             'read003\t147\tMULE_MUDRAVI2\t1093\t57\t7S28M65S\t=\t732\t-389\tTGATGACACGGTTTTGATTGTTTTTAAAAAGTAAAAACGTAATTGGTAACTACGATTGTATGACATGTTAATTTACTTAGAGTAAACGCATTAGATTGGA\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tNM:i:0\tMD:Z:28\tAS:i:28\tXS:i:20']

    answer = ['@read001\nATTATCATGTCATAAAATCGTAGTTACCAATTACGGTTTTC\n+\n12222222222222222222222222222222222222222',
              '@read003\nTCCAATCTAATGCGTTTACTCTAAGTAAATTAACATGTCATACAATCGTAGTTACCAATTACGTT\n+\n22222222222222222222222222222222222222222222222222222222222222222']

    assert preprocess.reverse_soft_clipped_tails_as_fastq(query) == answer
