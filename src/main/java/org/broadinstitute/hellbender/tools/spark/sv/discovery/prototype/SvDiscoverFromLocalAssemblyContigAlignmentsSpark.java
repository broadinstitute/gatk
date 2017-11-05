package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * This tool takes a SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads--currently the primary target),
 * searches for reads with--potentially complicated--split alignments indicating the presence of SV breakpoints.
 * In case of complex structural variations (cxSV), also outputs a custom file format containing how affected
 * reference segments are rearranged on the provided sample where the long reads are from.
 */
@CommandLineProgramProperties(summary="Parses a SAM file containing long reads alignments, and outputs cxSV rearrangements.",
        oneLineSummary="Parses a long read SAM file, and outputs cxSV rearrangements.",
        usageExample = "SvDiscoverFromLocalAssemblyContigAlignmentsSpark \\" +
                "-I /path/to/my/dir/localAssemblies.sam \\" +
                "-O /path/to/my/dir/outputDir \\" +
                "-R /path/to/my/reference/reference.2bit",
        omitFromCommandLine = true,
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class SvDiscoverFromLocalAssemblyContigAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(SvDiscoverFromLocalAssemblyContigAlignmentsSpark.class);

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
        final SAMFileHeader headerForReads = getHeaderForReads();
        final SAMSequenceDictionary refSequenceDictionary = headerForReads.getSequenceDictionary();

        final String sampleId = SVUtils.getSampleId(headerForReads);

        // filter alignments and split the gaps
        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed =
                FilterLongReadAlignmentsSAMSpark.filterByScore(reads, headerForReads, nonCanonicalChromosomeNamesFile, localLogger, 0.0)
                        .filter(lr -> lr.alignmentIntervals.size()>1).cache();

        // divert the long reads by their possible type of SV
        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes =
                divertReadsByPossiblyRawTypes(contigsWithAlignmentsReconstructed, localLogger);

        try {
            IOUtils.createDirectory(outputDir);
        } catch (final IOException x) {
            throw new UserException.CouldNotCreateOutputFile("Could not create file at path:" + outputDir + " due to " + x.getMessage(), x);
        }

        if (writeSAMFiles) {
            final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(headerForReads);
            contigsByPossibleRawTypes.forEach((k, v) -> writeSAM(v, k.name(), reads, headerBroadcast, outputDir, localLogger));
        }

        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(getReference());
        final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary = ctx.broadcast(refSequenceDictionary);

        dispatchJobs( contigsByPossibleRawTypes, referenceMultiSourceBroadcast, broadcastSequenceDictionary, sampleId);
    }

    private void dispatchJobs(final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes,
                              final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast,
                              final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary, final String sampleId) {

        new InsDelVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.InsDel), outputDir+"/"+RawTypes.InsDel.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        new SimpleStrandSwitchVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.Inv), outputDir+"/"+RawTypes.Inv.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        // here contigs supporting cpx having only 2 alignments could be handled easily together with ref block order switch ones
        final JavaRDD<AlignedContig> diffChrTrans =
                contigsByPossibleRawTypes.get(RawTypes.Cpx).filter(tig -> tig.alignmentIntervals.size() == 2);
        new SuspectedTransLocDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.DispersedDupOrMEI).union(diffChrTrans),
                        outputDir+"/"+ RawTypes.DispersedDupOrMEI.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);
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
                        contigWithOnlyOneConfigAnd2Aln.toString());
        return contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(0).referenceSpan.getContig()
                .equals(contigWithOnlyOneConfigAnd2Aln.alignmentIntervals.get(1).referenceSpan.getContig());
    }

    static boolean isLikelyInvBreakpointOrInsInv(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChr) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChr),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChr.toString());
        return contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(0).forwardStrand
                ^
                contigWithOnlyOneConfigAnd2AlnToSameChr.alignmentIntervals.get(1).forwardStrand;
    }

    private static boolean isSuggestingRefBlockOrderSwitch(final AlignedContig contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch) {
        Utils.validateArg(isSameChromosomeMapping(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.toString());
        Utils.validateArg( !isLikelyInvBreakpointOrInsInv(contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch),
                "assumption that input contig's 2 alignments map to the same chr is violated. \n" +
                        contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.toString());

        final AlignmentInterval one = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(0),
                                two = contigWithOnlyOneConfigAnd2AlnToSameChrWithoutStrandSwitch.alignmentIntervals.get(1);
        if (one.referenceSpan.contains(two.referenceSpan) || two.referenceSpan.contains(one.referenceSpan))
            return false;
        final List<AlignmentInterval> deOverlappedTempAlignments =
                ContigAlignmentsModifier.removeOverlap(one, two, null);

        return deOverlappedTempAlignments.get(0).referenceSpan.getStart() > deOverlappedTempAlignments.get(1).referenceSpan.getStart()
                == deOverlappedTempAlignments.get(0).forwardStrand;
    }

    private static boolean isLikelyCpx(final AlignedContig contigWithOnlyOneConfig) {
        Utils.validateArg(!contigWithOnlyOneConfig.hasEquallyGoodAlnConfigurations,
                "assumption that input contig has one unique best alignment configuration is violated: " +
                        contigWithOnlyOneConfig.toString());

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
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> cpxAndNot =
                RDDUtils.split(contigsWithOnlyOneBestConfig, SvDiscoverFromLocalAssemblyContigAlignmentsSpark::isLikelyCpx, true);

        contigsByRawTypes.put(RawTypes.Cpx, cpxAndNot._1);

        // long reads with only 1 best configuration and having only 2 alignments mapped to the same chromosome
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfigAnd2AIToSameChr = cpxAndNot._2;

        // divert away those with strand switch
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> strandSwitchAndNot =
                RDDUtils.split(contigsWithOnlyOneBestConfigAnd2AIToSameChr, SvDiscoverFromLocalAssemblyContigAlignmentsSpark::isLikelyInvBreakpointOrInsInv, true);

        contigsByRawTypes.put(RawTypes.Inv, strandSwitchAndNot._1);

        // 2 AI, same chr, no strand switch, then only 2 cases left
        final JavaRDD<AlignedContig> contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch = strandSwitchAndNot._2;

        // case 1: dispersed duplication, or MEI (that is, reference blocks seemingly switched their orders)
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> x =
                RDDUtils.split(contigsWithOnlyOneBestConfigAnd2AIToSameChrWithoutStrandSwitch, SvDiscoverFromLocalAssemblyContigAlignmentsSpark::isSuggestingRefBlockOrderSwitch, true);
        contigsByRawTypes.put(RawTypes.DispersedDupOrMEI, x._1);

        // case 2: no order switch: ins, del, or tandem dup
        contigsByRawTypes.put(RawTypes.InsDel, x._2);

        return contigsByRawTypes;
    }

    private static void writeSAM(final JavaRDD<AlignedContig> filteredContigs, final String rawTypeString,
                                 final JavaRDD<GATKRead> originalContigs, final Broadcast<SAMFileHeader> headerBroadcast,
                                 final String outputDir, final Logger toolLogger) {

        final Set<String> filteredReadNames = new HashSet<>( filteredContigs.map(tig -> tig.contigName).distinct().collect() );
        toolLogger.info(filteredReadNames.size() + " long reads indicating " + rawTypeString);
        final JavaRDD<SAMRecord> splitLongReads = originalContigs.filter(read -> filteredReadNames.contains(read.getName()))
                .map(read -> read.convertToSAMRecord(headerBroadcast.getValue()));
        SVFileUtils.writeSAMFile(outputDir+"/"+rawTypeString+".sam", splitLongReads.collect().iterator(), headerBroadcast.getValue(),
                false);
    }

}
