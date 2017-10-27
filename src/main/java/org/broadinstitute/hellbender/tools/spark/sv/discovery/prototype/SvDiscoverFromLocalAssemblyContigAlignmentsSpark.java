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
import org.broadinstitute.hellbender.tools.spark.sv.utils.FileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AssemblyContigAlignmentSignatureClassifier.RawTypes;

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
        final Broadcast<SAMFileHeader> headerBroadcast = ctx.broadcast(headerForReads);
        final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast = ctx.broadcast(getReference());
        final SAMSequenceDictionary refSequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary = ctx.broadcast(refSequenceDictionary);
        final String sampleId = SVUtils.getSampleId(headerForReads);

        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes =
                preprocess(reads, headerBroadcast, broadcastSequenceDictionary, nonCanonicalChromosomeNamesFile, outputDir, writeSAMFiles, localLogger);

        dispatchJobs(sampleId, contigsByPossibleRawTypes, referenceMultiSourceBroadcast, broadcastSequenceDictionary);
    }

    //==================================================================================================================

    private static EnumMap<RawTypes, JavaRDD<AlignedContig>> preprocess(final JavaRDD<GATKRead> reads,
                                                                        final Broadcast<SAMFileHeader> headerBroadcast,
                                                                        final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                        final String nonCanonicalChromosomeNamesFile,
                                                                        final String outputDir,
                                                                        final boolean writeSAMFiles,
                                                                        final Logger localLogger) {
        // filter alignments and split the gaps, hence the name "reconstructed"
        final JavaRDD<AlignedContig> contigsWithAlignmentsReconstructed =
                FilterLongReadAlignmentsSAMSpark.filterByScore(reads, headerBroadcast.getValue(), nonCanonicalChromosomeNamesFile, 0.0, localLogger)
                        .filter(lr -> lr.alignmentIntervals.size() > 1).cache();

        // classify and divert assembly contigs by their possible type of SV
        final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes =
                AssemblyContigAlignmentSignatureClassifier.classifyContigs(contigsWithAlignmentsReconstructed,
                        broadcastSequenceDictionary, localLogger);

        // prepare for output and optionally writes SAM files separately for contigs of different classes
        if ( !FileUtils.createDirToWriteTo(outputDir) )
            throw new UserException.CouldNotCreateOutputFile(outputDir,
                    "Could not create directory to output results", new IOException(""));

        if (writeSAMFiles) {
            contigsByPossibleRawTypes.forEach((k, v) -> writeSAM(v, k.name(), reads, headerBroadcast, outputDir, localLogger));
        }

        return contigsByPossibleRawTypes;
    }

    //==================================================================================================================

    private void dispatchJobs(final String sampleId,
                              final EnumMap<RawTypes, JavaRDD<AlignedContig>> contigsByPossibleRawTypes,
                              final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast,
                              final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary) {

        new InsDelVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.InsDel), outputDir+"/"+ RawTypes.InsDel.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        new SimpleStrandSwitchVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.IntraChrStrandSwitch), outputDir+"/"+ RawTypes.IntraChrStrandSwitch.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        new SuspectedTransLocDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.TandemDupOrMEIBkpt),
                        outputDir+"/"+ RawTypes.TandemDupOrMEIBkpt.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);
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
}
