from unittest import TestCase

from gcnvkernel.structs.metadata import SampleMetadataCollection
from gcnvkernel.postprocess.viterbi_segmentation import ViterbiSegmentationEngine
import gcnvkernel.io.io_metadata as io_metadata


class TestViterbiSegmentationEngine(TestCase):
    def test__viterbi_segments_generator(self):

        test_sub_dir = "/Users/gauthier/workspaces/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-postprocess/"
        ploidy_calls_path = test_sub_dir + "ploidy-calls/"
        model_shards = (test_sub_dir + "shard_0-model", test_sub_dir + "shard_1-model", test_sub_dir + "shard_2-model")
        calls_shards = (test_sub_dir + "shard_0-calls", test_sub_dir + "shard_1-calls", test_sub_dir + "shard_2-calls")
        sample_index = 0
        output_path = "."
        combined_intervals_vcf = test_sub_dir + "intervals.combined.vcf.gz"
        clustered_vcf = test_sub_dir + "../clustering/threeSamples.vcf.gz"

        sample_metadata_collection: SampleMetadataCollection = SampleMetadataCollection()
        io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
            sample_metadata_collection, ploidy_calls_path)

        viterbi_engine = ViterbiSegmentationEngine(
            model_shards, calls_shards, sample_metadata_collection, sample_index, output_path,
            combined_intervals_vcf, clustered_vcf)

        segments = list(viterbi_engine._viterbi_segments_generator())
        self.assertTrue(len(segments) == 11)
        ref_seg_0 = segments[0]
        self.assertTrue(ref_seg_0.num_points == 102)
        self.assertTrue(ref_seg_0.start == 68993)
        self.assertTrue(ref_seg_0.end == 986037)
        self.assertTrue(ref_seg_0.call_copy_number == 2)
        self.assertTrue(ref_seg_0.quality_some_called > ref_seg_0.quality_all_called)

        sample_index = 1
        viterbi_engine = ViterbiSegmentationEngine(
            model_shards, calls_shards, sample_metadata_collection, sample_index, output_path,
            combined_intervals_vcf, clustered_vcf)
        segments1 = list(viterbi_engine._viterbi_segments_generator())
        self.assertTrue(len(segments1) == len(segments))
        ref_seg_1 = segments1[0]
        self.assertTrue(ref_seg_1.num_points == ref_seg_0.num_points)
        self.assertTrue(ref_seg_1.start == ref_seg_0.start)
        self.assertTrue(ref_seg_1.end == ref_seg_0.end)
        self.assertTrue(ref_seg_1.call_copy_number == ref_seg_0.call_copy_number)
        self.assertTrue(ref_seg_1.quality_some_called != ref_seg_0.quality_some_called)


    def test_write_copy_number_segments(self):
        self.fail()

    def test__coalesce_seq_into_segments(self):
        self.fail()
