package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.FileUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.EnumMap;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;


/**
 * This tool takes a SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads--currently the primary target),
 * searches for reads with--potentially complicated--split alignments indicating the presence of SV breakpoints.
 * In case of complex structural variations (cxSV), also outputs a custom file format containing how affected
 * reference segments are rearranged on the provided sample where the long reads are from.
 */
@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs cxSV rearrangements.",
        oneLineSummary="Parses a long read SAM file, and outputs cxSV rearrangements.",
        usageExample = "InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark \\" +
                "-I /path/to/my/dir/localAssemblies.sam \\" +
                "-O /path/to/my/dir/outputDir \\" +
                "-R /path/to/my/reference/reference.2bit --fastaReference /path/to/my/reference/reference.fasta",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, human reference (hg19 or hg38) assumed when omitted",
            shortName = "nonCanoChrFile",
            fullName = "nonCanonicalChromosomeNamesFile", optional = true)
    public String nonCanonicalChromosomeNamesFile;

    @Argument(doc = "output directory", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Argument(doc = "output SAM files", shortName = "wSAM",
              fullName = "writeSAM", optional = true)
    private boolean writeSAMFiles = false;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();
        final SAMFileHeader header = getHeaderForReads();

        // filter alignments and split the gaps
        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed =
                InternalFilterLongReadAlignmentsSAMSpark.filterByScore(reads, header, nonCanonicalChromosomeNamesFile, localLogger)
                        .filter(lr -> lr.alignmentIntervals.size()>1).cache();

        // divert the long reads by their possible type of SV
        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes =
                divertReadsByPossiblyRawTypes(contigsWithAlignmentsReconstructed, localLogger);

        if ( !FileUtils.createDirToWriteTo(outputDir) )
            throw new GATKException("Could not create directory " + outputDir + " to write results to.");

        if (writeSAMFiles) {
            final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(header);
            contigsByPossibleRawTypes.forEach((k, v) -> writeSAM(v, k.name(), reads, headerBroadcast, outputDir, localLogger));
        }

        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(getReference());

        new InsDelVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.InsDel), outputDir+"/"+RawTypes.InsDel.name()+".vcf",
                        referenceMultiSourceBroadcast, discoverStageArgs.fastaReference, localLogger);

        new SimpleStrandSwitchVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.Inv), outputDir+"/"+RawTypes.Inv.name()+".vcf",
                        referenceMultiSourceBroadcast, discoverStageArgs.fastaReference, localLogger);
    }

    private enum RawTypes {
        Ambiguous, Inv, InsDel, DispersedDupOrMEI, Cpx;
    }

    private static boolean hasOnly2Alignments(final AlignedContig contigWithOnlyOneConfig) {
        return contigWithOnlyOneConfig.alignmentIntervals.size() == 2;
    }

    private static boolean isSameChromosomeMapping(final AlignedContig contigWithOnlyOneConfigAnd2Aln) {
        Utils.validateArg(hasOnly2Alignments(contigWithOnlyOneConfigAnd2Aln),
                "assumption that input contig has only 2 alignments is violated. \n" +
                        onErrorStringRepForAlignedContig(contigWithOnlyOneConfigAnd2Aln));
        return contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0).referenceSpan.getContig()
                .equals(contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1).referenceSpan.getContig());
    }

    static boolean isLikelyInvBreakpointOrInsInv(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChr) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChr),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        onErrorStringRepForAlignedContig(contigWithOnlyOneConfigAnd2AlnToSameChr));
        return contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0).forwardStrand
                ^
                contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1).forwardStrand;
    }

    private static boolean isSuggestingRefBlockOrderSwitch(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        onErrorStringRepForAlignedContig(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch));
        Utils.validateArg( !isLikelyInvBreakpointOrInsInv(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        onErrorStringRepForAlignedContig(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch));

        final AlignmentInterval intervalOne = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(0),
                                intervalTwo = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(1);

        return intervalOne.referenceSpan.getStart() > intervalTwo.referenceSpan.getStart() == intervalOne.forwardStrand;
    }

    private static boolean isLikelyCpx(final AlignedContig contigWithOnlyOneConfig) {
        Utils.validateArg(!contigWithOnlyOneConfig.hasEquallyGoodAlnConfigurations,
                "assumption that input contig has one unique best alignment configuration is violated: " +
                        onErrorStringRepForAlignedContig(contigWithOnlyOneConfig));

        return !hasOnly2Alignments(contigWithOnlyOneConfig) || !isSameChromosomeMapping(contigWithOnlyOneConfig);
    }

    private static EnumMap<RawTypes, JavaRDD<AlignedContig>> divertReadsByPossiblyRawTypes(final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed,
                                                                                           final Logger toolLogger) {

        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByRawTypes = new EnumMap<>(RawTypes.class);

        // long reads with more than 1 best configurations
        contigsByRawTypes.put(RawTypes.Ambiguous,
                contigsWithAlignmentsReconstructed.filter(tig -> tig.hasEquallyGoodAlnConfigurations));

        // long reads with only 1 best configuration
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfig =
                contigsWithAlignmentsReconstructed.filter(lr -> !lr.hasEquallyGoodAlnConfigurations).cache();

        // divert away those likely suggesting cpx sv (more than 2 alignments after gap split, or 2 alignments to diff chr)
        contigsByRawTypes.put(RawTypes.Cpx,
                contigsWithOnlyOneBestConfig.filter(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark::isLikelyCpx));

        // long reads with only 1 best configuration and having only 2 alignments mapped to the same chromosome
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfigAnd2AIToSameChr =
                contigsWithOnlyOneBestConfig
                        .filter(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark::hasOnly2Alignments)
                        .filter(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark::isSameChromosomeMapping).cache();

        // divert away those with strand switch
        contigsByRawTypes.put(RawTypes.Inv,
                contigsWithOnlyOneBestConfigAnd2AIToSameChr.filter(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark::isLikelyInvBreakpointOrInsInv));

        // 2 AI, same chr, no strand switch, then only 2 cases left
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch =
                contigsWithOnlyOneBestConfigAnd2AIToSameChr.filter(tig -> !isLikelyInvBreakpointOrInsInv(tig)).cache();

        // case 1: dispersed duplication, or MEI (that is, reference blocks seemingly switched their orders)
        contigsByRawTypes.put(RawTypes.DispersedDupOrMEI,
                contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch.filter(InternalSvDiscoverFromLocalAssemblyContigAlignmentsSpark::isSuggestingRefBlockOrderSwitch));

        // case 2: no order switch: ins, del, or tandem dup
        contigsByRawTypes.put(RawTypes.InsDel,
                contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch.filter(tig -> !isSuggestingRefBlockOrderSwitch(tig)));

        contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch.unpersist();
        contigsWithOnlyOneBestConfigAnd2AIToSameChr.unpersist();
        contigsWithOnlyOneBestConfig.unpersist();

        return contigsByRawTypes;
    }

    private static void writeSAM(final JavaRDD<AlignedContig> filteredContigs, final String rawTypeString,
                                 final JavaRDD<GATKRead> originalContigs, final Broadcast<SAMFileHeader> headerBroadcast,
                                 final String outputDir, final Logger toolLogger) {

        final Set<String> filteredReadNames = new HashSet<>( filteredContigs.map(tig -> tig.contigName).distinct().collect() );
        toolLogger.info(filteredReadNames.size() + " long reads indicating " + rawTypeString);
        final JavaRDD<SAMRecord> splitLongReads = originalContigs.filter(read -> filteredReadNames.contains(read.getName()))
                .map(read -> read.convertToSAMRecord(headerBroadcast.getValue()));
        FileUtils.writeSAMFile(splitLongReads.collect().iterator(), headerBroadcast.getValue(),
                outputDir+"/"+rawTypeString+".sam", false);
    }

    static String onErrorStringRepForAlignedContig(final AlignedContig contig) {
        return InternalFilterLongReadAlignmentsSAMSpark.formatContigInfo(
                new Tuple2<>(contig.contigName ,
                        contig.alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())));
    }
}
