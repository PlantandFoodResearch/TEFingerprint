#! /usr/bin/env python
import tempfile
import os
import pytest
import subprocess
from tefingerprint import preprocess


@pytest.mark.parametrize('query,answer',
                         [('AGATC', 'GATCT'),
                          ('aGtcn', 'NGACT')])
def test_reverse_complement(query, answer):
    assert preprocess.reverse_complement(query) == answer


@pytest.mark.parametrize('include_tails,minimum_tail,answer',
                         [(False, 38,
                           [{'element': 'MULE_1',
                             'name': 'read01_forward_dangler',
                             'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
                             'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
                            {'element': 'Gypsy_2',
                             'name': 'read02_forward_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
                             'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
                            {'element': 'MULE_3',
                             'name': 'read03_forward_dangler',
                             'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
                             'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
                            {'element': 'Gypsy_1',
                             'name': 'read04_forward_dangler',
                             'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
                             'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
                            {'element': 'MULE_2',
                             'name': 'read05_reverse_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
                             'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
                            {'element': 'Gypsy_3',
                             'name': 'read06_reverse_dangler',
                             'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
                             'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
                            {'element': 'MULE_1',
                             'name': 'read07_reverse_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
                             'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
                            {'element': 'Gypsy_2',
                             'name': 'read08_reverse_dangler',
                             'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
                             'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'}]),
                          (True, 38,
                           [{'element': 'MULE_1',
                             'name': 'read01_forward_dangler',
                             'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
                             'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
                            {'element': 'Gypsy_2',
                             'name': 'read02_forward_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
                             'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
                            {'element': 'MULE_3',
                             'name': 'read03_forward_dangler',
                             'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
                             'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
                            {'element': 'Gypsy_1',
                             'name': 'read04_forward_dangler',
                             'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
                             'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
                            {'element': 'MULE_2',
                             'name': 'read05_reverse_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
                             'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
                            {'element': 'Gypsy_3',
                             'name': 'read06_reverse_dangler',
                             'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
                             'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
                            {'element': 'MULE_1',
                             'name': 'read07_reverse_dangler',
                             'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
                             'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
                            {'element': 'Gypsy_2',
                             'name': 'read08_reverse_dangler',
                             'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
                             'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'},
                            {'element': 'MULE_3',
                             'name': 'read09_forward_tail',
                             'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ',
                             'sequence': 'TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT'},
                            {'element': 'Gypsy_3',
                             'name': 'read12_reverse_tail',
                             'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
                             'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'}]),
                          (True, 10,
                            [{'element': 'MULE_1',
                              'name': 'read01_forward_dangler',
                              'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
                              'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
                             {'element': 'Gypsy_2',
                              'name': 'read02_forward_dangler',
                              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
                              'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
                             {'element': 'MULE_3',
                              'name': 'read03_forward_dangler',
                              'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
                              'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
                             {'element': 'Gypsy_1',
                              'name': 'read04_forward_dangler',
                              'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
                              'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
                             {'element': 'MULE_2',
                              'name': 'read05_reverse_dangler',
                              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
                              'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
                             {'element': 'Gypsy_3',
                              'name': 'read06_reverse_dangler',
                              'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
                              'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
                             {'element': 'MULE_1',
                              'name': 'read07_reverse_dangler',
                              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
                              'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
                             {'element': 'Gypsy_2',
                              'name': 'read08_reverse_dangler',
                              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
                              'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'},
                             {'element': 'MULE_3',
                              'name': 'read09_forward_tail',
                              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ',
                              'sequence': 'TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT'},
                             {'element': 'MULE_2',
                              'name': 'read11_forward_short_tail',
                              'quality': 'HHHHHHHHHHH',
                              'sequence': 'AATTAAATGAA'},
                             {'element': 'Gypsy_3',
                              'name': 'read12_reverse_tail',
                              'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
                              'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'},
                             {'element': 'MULE_2',
                              'name': 'read14_reverse_short_tail',
                              'quality': '@CCFFFDDHF',
                              'sequence': 'TCCACTGTTT'}])
                          ])
def test_extract_informative_reads(include_tails, minimum_tail, answer):
    """
    Tests read extraction parsing and filtering.
    """
    data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/testPreprocessInput-2017-07-28.bam'
    query = list(preprocess.extract_informative_reads(data_path, include_tails=include_tails, minimum_tail=minimum_tail))
    assert query == answer


def test_capture_elements_and_write_fastq():
    query = [{'element': 'MULE_1',
              'name': 'read01_forward_dangler',
              'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
              'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
             {'element': 'Gypsy_2',
              'name': 'read02_forward_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
              'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
             {'element': 'MULE_3',
              'name': 'read03_forward_dangler',
              'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
              'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
             {'element': 'Gypsy_1',
              'name': 'read04_forward_dangler',
              'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
              'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
             {'element': 'MULE_2',
              'name': 'read05_reverse_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
              'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
             {'element': 'Gypsy_3',
              'name': 'read06_reverse_dangler',
              'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
              'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
             {'element': 'MULE_1',
              'name': 'read07_reverse_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
              'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
             {'element': 'Gypsy_2',
              'name': 'read08_reverse_dangler',
              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
              'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'},
             {'element': 'MULE_3',
              'name': 'read09_forward_tail',
              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ',
              'sequence': 'TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT'},
             {'element': 'Gypsy_3',
              'name': 'read12_reverse_tail',
              'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
              'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'}]

    answer_elements = {'read01_forward_dangler': 'MULE_1',
                       'read02_forward_dangler': 'Gypsy_2',
                       'read03_forward_dangler': 'MULE_3',
                       'read04_forward_dangler': 'Gypsy_1',
                       'read05_reverse_dangler': 'MULE_2',
                       'read06_reverse_dangler': 'Gypsy_3',
                       'read07_reverse_dangler': 'MULE_1',
                       'read08_reverse_dangler': 'Gypsy_2',
                       'read09_forward_tail': 'MULE_3',
                       'read12_reverse_tail': 'Gypsy_3'}

    answer_fastq = """@read01_forward_dangler
TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT
+
HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE
@read02_forward_dangler
TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH
@read03_forward_dangler
GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT
+
CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD
@read04_forward_dangler
ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG
+
@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B
@read05_reverse_dangler
TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE
@read06_reverse_dangler
TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA
+
@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE
@read07_reverse_dangler
CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=
@read08_reverse_dangler
ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA
+
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF
@read09_forward_tail
TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT
+
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ
@read12_reverse_tail
TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT
+
@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ
"""

    with tempfile.NamedTemporaryFile() as fastq:
        elements = preprocess.capture_elements_and_write_fastq(query, fastq.name)
        assert elements == answer_elements

        with open(fastq.name, 'r') as handle:
            assert handle.read() == answer_fastq
    fastq.close()


def test_tag_danglers():
    pass
