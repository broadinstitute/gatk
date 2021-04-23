from unittest import TestCase
import os

from gcnvkernel.structs.metadata import SampleMetadataCollection
from gcnvkernel.postprocess.viterbi_segmentation import ViterbiSegmentationEngine
import gcnvkernel.io.io_metadata as io_metadata


class test_viterbiSegmentationEngine(TestCase):
    def test__viterbi_segments_generator(self):

        current_dir = os.getcwd()
        #for GATK PythonUnitTestRunner/Java tests
        test_sub_dir = current_dir + "/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        # for running in IntelliJ/Python tests
        # test_sub_dir = current_dir + "/../../../../../../../../src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        ploidy_calls_path = test_sub_dir + "ploidy-calls/"
        model_shards = (test_sub_dir + "shard_0-model", test_sub_dir + "shard_1-model", test_sub_dir + "shard_2-model")
        calls_shards = (test_sub_dir + "shard_0-calls", test_sub_dir + "shard_1-calls", test_sub_dir + "shard_2-calls")
        sample_index = 0
        intervals_vcf = test_sub_dir + "intervals_output_SAMPLE_000.vcf.gz"
        output_path = "."
        clustered_vcf = test_sub_dir + "../clustering/threeSamples.vcf.gz"

        sample_metadata_collection: SampleMetadataCollection = SampleMetadataCollection()
        io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
            sample_metadata_collection, ploidy_calls_path)

        viterbi_engine = ViterbiSegmentationEngine(
            model_shards, calls_shards, sample_metadata_collection, sample_index, output_path,
            intervals_vcf, clustered_vcf)

        segments = list(viterbi_engine._viterbi_segments_generator())
        self.assertTrue(len(segments) == 6) #to match the number of VCs in the clustering VCF
        ref_seg_0 = segments[0]
        self.assertTrue(ref_seg_0.num_points == 1)
        self.assertTrue(ref_seg_0.start == 230925)
        self.assertTrue(ref_seg_0.end == 231288)
        self.assertTrue(ref_seg_0.call_copy_number == 2)
        self.assertTrue(ref_seg_0.quality_some_called - ref_seg_0.quality_all_called < 0.01) #should be the same, but numerical precision

        sample_index = 1
        intervals_vcf = test_sub_dir + "intervals_output_SAMPLE_001.vcf.gz"
        viterbi_engine = ViterbiSegmentationEngine(
            model_shards, calls_shards, sample_metadata_collection, sample_index, output_path,
            intervals_vcf, clustered_vcf)
        segments1 = list(viterbi_engine._viterbi_segments_generator())
        self.assertTrue(len(segments1) == len(segments))
        del_seg = segments1[0]
        self.assertTrue(del_seg.num_points == ref_seg_0.num_points)
        self.assertTrue(del_seg.start == ref_seg_0.start)
        self.assertTrue(del_seg.end == ref_seg_0.end)
        self.assertTrue(del_seg.call_copy_number == 0)
        self.assertTrue(del_seg.quality_some_called - del_seg.quality_all_called < 0.01)

    def test__sample_name_mismatch(self):
        current_dir = os.getcwd()
        test_sub_dir = current_dir + "/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        ploidy_calls_path = test_sub_dir + "ploidy-calls/"
        model_shards = (test_sub_dir + "shard_0-model", test_sub_dir + "shard_1-model", test_sub_dir + "shard_2-model")
        calls_shards = (test_sub_dir + "shard_0-calls", test_sub_dir + "shard_1-calls", test_sub_dir + "shard_2-calls")
        sample_index = 0
        output_path = "."
        combined_intervals_vcf = test_sub_dir + "intervals_output_SAMPLE_000.vcf.gz"
        # any test VCF that doesn't have SAMPLE_000
        clustered_vcf = test_sub_dir + "../../walkers/sv/JointGermlineCNVSegmentation/NA20533.fragmented.segments.vcf.gz"

        sample_metadata_collection: SampleMetadataCollection = SampleMetadataCollection()
        io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
             sample_metadata_collection, ploidy_calls_path)

        with self.assertRaises(AssertionError):
            viterbi_engine = ViterbiSegmentationEngine(
                model_shards, calls_shards, sample_metadata_collection, sample_index, output_path,
                combined_intervals_vcf, clustered_vcf)