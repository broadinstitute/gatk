from unittest import TestCase
import os
import gcnvkernel.io.io_vcf_parsing as io


class test_io_vcf_parsing(TestCase):
    def test_read_sample_segments_and_calls(self):

        current_dir = os.getcwd()
        #for GATK PythonUnitTestRunner/Java tests
        test_sub_dir = current_dir + "/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        # for running in IntelliJ/Python tests
        # test_sub_dir = current_dir + "/../../../../../../../../src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"

        clustered_vcf = test_sub_dir + 'clustered.1000G.vcf.gz'
        pesky_intervals_vcf = test_sub_dir + 'genotyped-intervals-HG00099.mapped.ILLUMINA.bwa.GBR.exome.20130415.bam.cram.vcf.gz'
        pesky_sample_name = 'HG00099'
        contig = '14'
        debug_path = io.read_sample_segments_and_calls(pesky_intervals_vcf, clustered_vcf, pesky_sample_name, contig)
        self.assertTrue(len(debug_path) == 12)  #should match number of chr 14 lines in clustered VCF (12)

        clustered_vcf = current_dir + '/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/clustering/threeSamples.vcf.gz'
        intervals_vcf = test_sub_dir + 'intervals_output_SAMPLE_000.vcf.gz'
        sample_name = 'SAMPLE_000'

        contig = "1"
        path = io.read_sample_segments_and_calls(intervals_vcf, clustered_vcf, sample_name, contig)
        # no segments on contig 1
        self.assertTrue(len(path) == 0)

        contig = "2"
        path = io.read_sample_segments_and_calls(intervals_vcf, clustered_vcf, sample_name, contig)
        self.assertTrue(len(path) == 2)

        sample1_intervals_vcf = test_sub_dir + 'intervals_output_SAMPLE_001.vcf.gz'
        sample1_name = 'SAMPLE_001'
        contig = "2"
        path1 = io.read_sample_segments_and_calls(sample1_intervals_vcf, clustered_vcf, sample1_name, contig)
        # all samples should have the same number of intervals
        self.assertTrue(len(path) == len(path1))
