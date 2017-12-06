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
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AssemblyContigAlignmentSignatureClassifier.RawTypes;

/**
 * Parse aligned contigs and call structural variants.
 *
 * <p>This tool takes a file containing the alignments of assembled contigs and searches for reads with
 * split alignments indicating the presence of SV breakpoints. The alignment signatures of the split alignments are
 * analyzed to determine the type of structural variation and written to a VCF file.</p>
 * <p>This is a prototyping/debugging tool and it is probably not generally useful to most users.</p>
 * <p>The input file is typically the SAM file produced by FindBreakpointEvidenceSpark.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of assembled contigs or long reads aligned to reference.</li>
 *     <li>The reference to which the contigs have been aligned.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>Text files describing the discovered structural variants and complex structural variants in the specified output directory.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk SvDiscoverFromLocalAssemblyContigAlignmentsSpark \
 *     -I assemblies.sam \
 *     -R reference.2bit \
 *     -O output_directory
 * </pre>
 *
 * <h3>Notes</h3>
 * <p>The reference is broadcast by Spark, and must therefore be a 2bit file due to current restrictions.</p>
 */
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Parse aligned contigs and call complex structural variants.",
        summary =
        "This tool takes a file containing the alignments of assembled contigs and searches for reads with\n" +
        " split alignments indicating the presence of SV breakpoints. The alignment signatures of the split alignments are" +
        " analyzed to determine the type of structural variation and written to a VCF file.",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class SvDiscoverFromLocalAssemblyContigAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(SvDiscoverFromLocalAssemblyContigAlignmentsSpark.class);

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection
            discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, human reference (hg19 or hg38) assumed when omitted",
            shortName = "alt-tigs",
            fullName = "non-canonical-contig-names-file", optional = true)
    public String nonCanonicalChromosomeNamesFile;

    @Argument(doc = "output directory for outputting VCF files for each type of variant, and if enabled, the signaling assembly contig's alignments",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Argument(doc = "output SAM files", fullName = "write-sam", optional = true)
    private boolean writeSAMFiles;

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

        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                preprocess(reads, headerBroadcast, broadcastSequenceDictionary, nonCanonicalChromosomeNamesFile, outputDir, writeSAMFiles, localLogger);

        dispatchJobs(sampleId, outputDir, contigsByPossibleRawTypes, referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger);
    }

    //==================================================================================================================

    public static EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> preprocess(final JavaRDD<GATKRead> reads,
                                                                                               final Broadcast<SAMFileHeader> headerBroadcast,
                                                                                               final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                                               final String nonCanonicalChromosomeNamesFile,
                                                                                               final String outputDir,
                                                                                               final boolean writeSAMFiles,
                                                                                               final Logger localLogger) {
        // filter alignments and split the gaps, hence the name "reconstructed"
        final JavaRDD<AlignedContig> contigsWithChimericAlignmentsReconstructed =
                FilterLongReadAlignmentsSAMSpark
                        .filterByScore(reads, headerBroadcast.getValue(), nonCanonicalChromosomeNamesFile, 0.0, localLogger)
                        .filter(lr -> lr.alignmentIntervals.size() > 1).cache();
        localLogger.info( contigsWithChimericAlignmentsReconstructed.count() +
                " contigs with chimeric alignments potentially giving SV signals.");

        // classify assembly contigs by their possible type of SV based on studying alignment signature
        final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes =
                AssemblyContigAlignmentSignatureClassifier.classifyContigs(contigsWithChimericAlignmentsReconstructed, broadcastSequenceDictionary, localLogger);

        debug(reads, contigsByPossibleRawTypes, headerBroadcast, outputDir, writeSAMFiles, localLogger);

        return contigsByPossibleRawTypes;
    }

    //==================================================================================================================

    // TODO: 11/21/17 insertion mappings are dropped here in this implementation, must get them back for ticket #3647

    /**
     * Sends assembly contigs classified based on their alignment signature to
     * a corresponding breakpoint location inference unit.
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#InsDel},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#IntraChrStrandSwitch},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#MappedInsertionBkpt}, and
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Cpx} (PR to be reviewed and merged)
     * each will have its own VCF output in the directory specified in {@link #outputDir},
     * whereas
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Incomplete},
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#Ambiguous}, and
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes#MisAssemblySuspect}
     * currently DO NOT generate any VCF yet.
     * However, if flag {@link #writeSAMFiles} is turned on, alignments of all contigs that are classified to be any of
     * {@link AssemblyContigAlignmentSignatureClassifier.RawTypes}
     * will be extracted and put in SAM files in {@link #outputDir} too.
     */
    public static void dispatchJobs(final String sampleId, final String outputDir,
                                    final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes,
                                    final Broadcast<ReferenceMultiSource> referenceMultiSourceBroadcast,
                                    final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                    final Logger localLogger) {

        new InsDelVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.InsDel).map(decoratedTig -> decoratedTig.contig),
                        outputDir +"/"+ RawTypes.InsDel.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        new SimpleStrandSwitchVariantDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.IntraChrStrandSwitch).map(decoratedTig -> decoratedTig.contig),
                        outputDir +"/"+ RawTypes.IntraChrStrandSwitch.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);

        new SuspectedTransLocDetector()
                .inferSvAndWriteVCF(contigsByPossibleRawTypes.get(RawTypes.MappedInsertionBkpt).map(decoratedTig -> decoratedTig.contig),
                        outputDir +"/"+ RawTypes.MappedInsertionBkpt.name()+".vcf",
                        referenceMultiSourceBroadcast, broadcastSequenceDictionary, localLogger, sampleId);
    }

    //==================================================================================================================
    
    private static void debug(final JavaRDD<GATKRead> reads, 
                              final EnumMap<RawTypes, JavaRDD<AssemblyContigWithFineTunedAlignments>> contigsByPossibleRawTypes, 
                              final Broadcast<SAMFileHeader> headerBroadcast,
                              final String outputDir, final boolean writeSAMFiles, final Logger localLogger) {
        try {
            IOUtils.createDirectory(outputDir);
            if (writeSAMFiles) {
                contigsByPossibleRawTypes.forEach((k, v) -> writeSAM(v, k.name(), reads, headerBroadcast, outputDir, localLogger));
            }
        } catch (final IOException x) {
            throw new UserException.CouldNotCreateOutputFile("Could not create file at path:" + outputDir + 
                    " due to " + x.getMessage(), x);
        }
    }
    
    private static void writeSAM(final JavaRDD<AssemblyContigWithFineTunedAlignments> filteredContigs, final String rawTypeString,
                                 final JavaRDD<GATKRead> originalContigs, final Broadcast<SAMFileHeader> headerBroadcast,
                                 final String outputDir, final Logger toolLogger) {

        final Set<String> filteredReadNames = new HashSet<>( filteredContigs.map(decoratedTig -> decoratedTig.contig).map(tig -> tig.contigName).distinct().collect() );
        toolLogger.info(filteredReadNames.size() + " contigs indicating " + rawTypeString);
        final JavaRDD<SAMRecord> splitLongReads = originalContigs.filter(read -> filteredReadNames.contains(read.getName()))
                .map(read -> read.convertToSAMRecord(headerBroadcast.getValue()));
        SVFileUtils.writeSAMFile(outputDir+"/"+rawTypeString+".sam", splitLongReads.collect().iterator(), headerBroadcast.getValue(),
                false);
    }
}
