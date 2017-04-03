package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.bwa.BwaSparkEngine;
import org.broadinstitute.hellbender.tools.spark.sv.ContainsKmerReadFilterSpark;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.ReadFilterSparkifier;
import org.broadinstitute.hellbender.tools.spark.utils.ReadTransformerSparkifier;
import org.broadinstitute.hellbender.transformers.BaseQualityClipReadTransformer;
import org.broadinstitute.hellbender.transformers.BaseQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.DUSTReadTransformer;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

/**
 * This Spark tool is the first step in the PathSeq sample processing pipeline. It takes in a BAM file
 * and filters the reads based on their sequences, base qualities, flags, and attributes. It then removes reads that
 * are sufficiently similar to the provided host organism (e.g. human) reference sequence.
 *
 * Filtering steps:
 *  1) Remove secondary, supplementary, and failed vendor quality check reads
 *  2) Remove optical duplicates
 *  3) Mask repetitive sequences with 'N' and base quality --dustPhred using symmetric DUST
 *  4) Hard clip read ends using base qualities
 *  5) Remove reads shorter than --minClippedReadLength
 *  6) Mask bases whose Phred score is less than --minBaseQuality with 'N'
 *  7) Remove reads whose fraction of bases that are 'N' is greater than --maxAmbiguousBaseFraction
 *  8) Remove reads containing one or more kmers from --kmerLibraryPath
 *  9) Remove reads that align to the host reference
 *
 * The tool assumes the BAM file is unaligned but will still work on aligned BAM files. However, it will ignore
 * any previous alignment information and therefore may repeat computations by re-aligning to the host reference.
 * All attributes except for RG will are discarded.
 *
 * The user must supply an indexed host FASTA reference and the host kmer library generated using PathSeqKmerSpark.
 *
 * The output is a BAM containing non-host reads ready to be aligned to a pathogen reference.
 */
@CommandLineProgramProperties(summary = "Read preprocessing and host organism filtering on reads from a BAM file",
        oneLineSummary = "PathSeqFilter on Spark",
        programGroup = SparkProgramGroup.class)
public final class PathSeqFilterSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;
    private SAMFileHeader header;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false)
    public String OUTPUT_PATH;

    @Argument(doc = "Path to kmer library generated with PathSeqKmerSpark",
            fullName="kmerLibraryPath",
            shortName="kLibPath",
            optional=false)
    public String KMER_LIB_PATH;

    @Argument(doc = "Path to indexed host reference FASTA",
            fullName="hostReference",
            shortName="hRef",
            optional=false)
    public String HOST_REF_PATH;

    @Argument(doc = "Keep only clipped reads with length at least equal to the specified value",
            fullName = "minClippedReadLength",
            shortName = "minClipLen",
            optional=true)
    public int MIN_READ_LENGTH = 31;

    @Argument(doc = "Max allowable fraction of ambiguous bases",
            fullName="maxAmbiguousBaseFraction",
            shortName="maxAmbigFrac",
            optional=true)
    public float FRAC_N_THRESHOLD = 0.05f;

    @Argument(doc = "Bases below this read quality will be transformed into 'N's",
            fullName="minBaseQuality",
            shortName="minBaseQual",
            optional=true)
    public int QUAL_PHRED_THRESH = 15;

    @Argument(doc = "Size of kmers in the kmer library",
            fullName="kmerSize",
            shortName="kSize",
            optional=true)
    public int KMER_SIZE = 31;

    @Argument(doc = "Quality score trimmer threshold",
            fullName = "readTrimmerThreshold",
            shortName = "readTrimThresh",
            optional = true)
    public int READ_TRIM_THRESH = 15;

    @Argument(doc = "Base quality to assign DUST masked bases",
            fullName="dustPhred",
            shortName="dustP",
            optional=true)
    public short DUST_MASK = 15;

    @Argument(doc = "DUST window size",
            fullName="dustWindowSize",
            shortName="dustW",
            optional=true)
    public int DUST_W = 64;

    @Argument(doc = "DUST score threshold",
            fullName="dustTScore",
            shortName="dustT",
            optional=true)
    public float DUST_T = 20.0f;

    @Argument(doc = "the bwa mem index image file name that you've distributed to each executor",
            fullName = "bwamemIndexImage")
    private String indexImageFile;

    @ArgumentCollection
    private OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    public boolean requiresReads() { return true; }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();

        //Filter secondary/supplementary reads and reads that fail the vendor quality check
        final JavaRDD<GATKRead> primaryReads = reads.filter(read -> !(read.isSecondaryAlignment() || read.failsVendorQualityCheck() || read.isSupplementaryAlignment()));
        logger.info("Loaded " + reads.count() + " reads");

        //Mark and filter optical duplicates
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null);
        final JavaRDD<GATKRead> markedReads = MarkDuplicatesSpark.mark(primaryReads, getHeaderForReads(),MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES, finder, getRecommendedNumReducers());
        final JavaRDD<GATKRead> markedFilteredReads = markedReads.filter(new ReadFilterSparkifier(new MarkedOpticalDuplicateReadFilter()));
        logger.info("Reads remaining after de-duplication: " + markedFilteredReads.count());

        //Apply DUST masking
        final JavaRDD<GATKRead> readsDUSTMasked = markedFilteredReads.map(new ReadTransformerSparkifier(new DUSTReadTransformer(DUST_MASK,DUST_W,DUST_T)));

        //Apply base quality hard clipping
        final JavaRDD<GATKRead> readsClipped = readsDUSTMasked.map(new ReadTransformerSparkifier(new BaseQualityClipReadTransformer(READ_TRIM_THRESH)));

        //Filter reads with less than MIN_READ_LENGTH bases
        final JavaRDD<GATKRead> readsLengthFiltered = readsClipped.filter(new ReadFilterSparkifier(new ReadLengthReadFilter(MIN_READ_LENGTH,Integer.MAX_VALUE)));
        logger.info("Reads remaining after clipping: " + readsLengthFiltered.count());

        //Change low-quality bases to 'N'
        final JavaRDD<GATKRead> readsBQFiltered = readsLengthFiltered.map(new ReadTransformerSparkifier(new BaseQualityReadTransformer(QUAL_PHRED_THRESH)));

        //Filter reads with too many 'N's
        final JavaRDD<GATKRead> readsAmbigFiltered = readsBQFiltered.filter(new ReadFilterSparkifier(new AmbiguousBaseReadFilter(FRAC_N_THRESHOLD)));
        logger.info("Reads remaining after ambiguous base filtering: " + readsAmbigFiltered.count());

        //Load Kmer hopscotch set and filter reads containing > 0 matching kmers
        final JavaRDD<GATKRead> readsKmerFiltered = doKmerFiltering(ctx,readsAmbigFiltered);
        logger.info("Reads remaining after kmer filtering: " + readsKmerFiltered.count());

        //Filter unpaired reads
        final JavaRDD<GATKRead> readsFilteredPaired = retainPairs(readsKmerFiltered);
        logger.info("Reads remaining after unpaired filtering: " + readsFilteredPaired.count());

        //BWA filtering against user-specified host organism reference
        header = getHeaderForReads();
        final JavaRDD<GATKRead> readsAligned = doHostBWA(ctx, header, readsFilteredPaired);

        //Get unmapped reads (note these always come in pairs)
        //TODO: retain read pairs by alignment score instead of flags
        final JavaRDD<GATKRead> readsNonHost = readsAligned.filter(read -> read.isUnmapped() && read.mateIsUnmapped());
        logger.info("Reads remaining after BWA filtering: " + readsFilteredPaired.count());

        //TODO: repeat BWA with seed size 11

        //Write output
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        try {
            ReadsSparkSink.writeReads(ctx, OUTPUT_PATH, null, readsNonHost, header, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (final IOException e) {
            throw new GATKException("Unable to write bam",e);
        }
    }

    private JavaRDD<GATKRead> doHostBWA(final JavaSparkContext ctx, final SAMFileHeader readsHeader, final JavaRDD<GATKRead> reads) {

        final BwaSparkEngine engine = new BwaSparkEngine(ctx, indexImageFile, getHeaderForReads(), getReferenceSequenceDictionary());
        final GCSOptions gcsOptions = getAuthenticatedGCSOptions(); // null if we have no api key
        final ReferenceMultiSource hostReference = new ReferenceMultiSource(gcsOptions, HOST_REF_PATH, getReferenceWindowFunction());
        final SAMSequenceDictionary hostRefDict = hostReference.getReferenceSequenceDictionary(header.getSequenceDictionary());
        readsHeader.setSequenceDictionary(hostRefDict);
        return engine.align(reads);
    }

    @SuppressWarnings("unchecked")
    private JavaRDD<GATKRead> doKmerFiltering(final JavaSparkContext ctx, final JavaRDD<GATKRead> reads) {

        final PipelineOptions options = getAuthenticatedGCSOptions();
        Input input = new Input(BucketUtils.openFile(KMER_LIB_PATH));
        Kryo kryo=new Kryo();
        kryo.setReferences(false);

        Set<SVKmer> kmerLibSet = (HopscotchSet<SVKmer>)kryo.readClassAndObject(input);

        return reads.filter(new ContainsKmerReadFilterSpark(ctx.broadcast(kmerLibSet),KMER_SIZE));
    }

    private static JavaRDD<GATKRead> retainPairs(JavaRDD<GATKRead> reads) {
        JavaRDD<GATKRead> pairedReads = reads.groupBy(read -> read.getName())
                                            .values()
                                            .filter(p -> {
                                                final Iterator<GATKRead> itr = p.iterator();
                                                itr.next();
                                                return itr.hasNext();})
                                            .flatMap(pairedRead -> pairedRead.iterator());
        return pairedReads;
    }

    /**
     * This version of retainPairs() groups read name Strings instead of the GATKReads themselves for faster shuffling.
     */
    private static JavaRDD<GATKRead> retainPairs_ShuffleByName(final JavaSparkContext ctx, JavaRDD<GATKRead> reads) {

        //Persist original reads since we're about to perform operations on it but need the original data later
        reads.cache();

        Collection<String> unpairedReadNames = reads.map(GATKRead::getName)
                .groupBy(read -> read)
                .values()
                .filter(p -> {
                    final Iterator<String> itr = p.iterator();
                    itr.next();
                    return !itr.hasNext();})
                .flatMap(unpairedReadName -> unpairedReadName.iterator())
                .collect();

        Broadcast<Collection<String>> unpairedReadNamesBroadcast = ctx.broadcast(unpairedReadNames);
        //The following line is quite slow when run locally
        return reads.filter(read -> !unpairedReadNamesBroadcast.value().contains(read.getName()));
    }

}
