import logging
import argparse
import gcnvkernel

logger = logging.getLogger("segment_gcnv_calls")
gcnvkernel.cli_commons.set_logging_config()

parser = argparse.ArgumentParser(description="gCNV segmentation tool",
                                 formatter_class=gcnvkernel.cli_commons.GCNVHelpFormatter)

# set theano flags and reload theano
gcnvkernel.cli_commons.set_theano_flags(parser)

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--ploidy_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="The path to the results of ploidy determination tool")

group.add_argument("--model_shards",
                   type=str,
                   required=True,
                   nargs='+',  # one or more
                   default=argparse.SUPPRESS,
                   help="List of coverage model shards (in genomic sequence order)")

group.add_argument("--calls_shards",
                   type=str,
                   required=True,
                   nargs='+',  # one or more
                   default=argparse.SUPPRESS,
                   help="List of gCNV calls shards (in genomic sequence order)")

group.add_argument("--output_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write segmented calls")

group.add_argument("--sample_index",
                   type=int,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Sample index to process")

if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()

    # load read depth and ploidy metadata
    logger.info("Loading ploidy calls...")
    logger.debug("Ploidy calls path: {0}".format(args.ploidy_calls_path))
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    gcnvkernel.io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection, args.ploidy_calls_path)

    # instantiate Viterbi segmentation engine
    logger.info("Instantiating the Viterbi segmentation engine...")
    logger.debug("Model shards path(s): {0}".format(repr(args.model_shards)))
    logger.debug("Calls shards path(s): {0}".format(repr(args.model_shards)))

    viterbi_engine = gcnvkernel.ViterbiSegmentationEngine(
        args.model_shards, args.calls_shards, sample_metadata_collection, args.sample_index, args.output_path)
    viterbi_engine.write_copy_number_segments()
