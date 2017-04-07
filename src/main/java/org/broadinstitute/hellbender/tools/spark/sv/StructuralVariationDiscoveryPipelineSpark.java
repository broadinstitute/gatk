package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.collections4.IterableUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tool to run the sv pipeline up to and including variant discovery
 */
@CommandLineProgramProperties(summary="Master tool to run the structural variation discovery pipeline",
        oneLineSummary="Master tool to run the structural variation discovery pipeline",
        programGroup = StructuralVariationSparkProgramGroup.class)
public class StructuralVariationDiscoveryPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(StructuralVariationDiscoveryPipelineSpark.class);


    @Argument(doc = "sam file for aligned contigs", shortName = "contigSAMFile",
            fullName = "contigSAMFile")
    private String outputSAM;

    @Argument(doc = "filename for output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection reAlignmentStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();


    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    public boolean requiresReference() {return true;}

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final SAMFileHeader header = getHeaderForReads();
        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();

        // gather evidence, run assembly, and align
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList =
                FindBreakpointEvidenceSpark
                        .gatherEvidenceAndWriteContigSamFile(ctx, pipelineOptions, reAlignmentStageArgs, header, getUnfilteredReads(),
                                outputSAM, localLogger);
        if (alignedAssemblyOrExcuseList.isEmpty()) return;

        // parse the contig alignments and extract necessary information
        @SuppressWarnings("unchecked")
        final JavaRDD<AlignedAssembly.AlignedContig> parsedAlignments = new InMemoryAlignmentParser().loadAndFilterAndParseContigAlignments(ctx, localLogger, header, alignedAssemblyOrExcuseList);
        if(parsedAlignments.isEmpty()) return;

        // discover variants and write to vcf
        DiscoverVariantsFromContigAlignmentsSAMSpark
                .discoverVariantsAndWriteVCF(parsedAlignments, discoverStageArgs.fastaReference,
                        ctx.broadcast(getReference()), pipelineOptions, vcfOutputFileName, localLogger);
    }

    public static final class InMemoryAlignmentParser implements AssemblyAlignmentParser {
        private static final long serialVersionUID = 1L;

        @SuppressWarnings("unchecked")
        public JavaRDD<AlignedAssembly.AlignedContig> loadAndFilterAndParseContigAlignments(final JavaSparkContext ctx,
                                                                                            final Logger toolLogger,
                                                                                            Object... localArguments) {

            final SAMFileHeader header = (SAMFileHeader) localArguments[0];
            final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList = (List<AlignedAssemblyOrExcuse>) localArguments[1];

            // here we have two options, one is going through the route "BwaMemAlignment -> SAM -> GATKRead -> SAM -> AlignmentInterval"
            //                           which is the route if the discovery pipeline is run by "FindBreakpointEvidenceSpark -> write sam file -> load sam file -> DiscoverVariantsFromContigAlignmentsSAMSpark"
            //                         , the other is to go directly "BwaMemAlignment -> AlignmentInterval", which is faster but not used here.
            return ctx.parallelize( viaSAMRoute(alignedAssemblyOrExcuseList, header, toolLogger) );
        }

        @VisibleForTesting
        public List<AlignedAssembly.AlignedContig> viaSAMRoute(final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                                               final SAMFileHeader header, final Logger toolLogger) {

            final SAMFileHeader cleanHeader = new SAMFileHeader(header.getSequenceDictionary());
            final List<AlignedAssembly.AlignedContig> parsedAlignments
                    = AlignedAssemblyOrExcuse.turnIntoSAMRecordsForEachAssembly(cleanHeader, alignedAssemblyOrExcuseList).stream()
                    .filter(samRecords -> Utils.stream(samRecords).anyMatch(sam -> !sam.getReadUnmappedFlag()))
                    .map(iterable -> Utils.stream(iterable).filter(sam -> !sam.getNotPrimaryAlignmentFlag()).collect(Collectors.toList()))
                    .filter(iterable -> !IterableUtils.isEmpty(iterable))
                    .map(iterable -> DiscoverVariantsFromContigAlignmentsSAMSpark.SAMFormattedContigAlignmentParser.parseReadsAndBreakGaps(Utils.stream(iterable).map(SAMRecordToGATKReadAdapter::new).collect(Collectors.toList()), header, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, toolLogger))
                    .collect(Collectors.toList());

            return parsedAlignments.isEmpty() ? Collections.emptyList() : parsedAlignments;
        }

        /**
         * Filter out "failed" assemblies and turn the alignments of contigs into custom {@link AlignedAssembly.AlignmentInterval} format.
         */
        @VisibleForTesting
        public static List<AlignedAssembly.AlignedContig> directToAlignmentInterval(final Iterable<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseIterable,
                                                                                    final List<String> refNames, final SAMFileHeader header) {

            return Utils.stream(alignedAssemblyOrExcuseIterable)
                    .filter(alignedAssemblyOrExcuse -> alignedAssemblyOrExcuse.getErrorMessage() == null)
                    .map(alignedAssembly -> forEachAssemblyNotExcuse(alignedAssembly, refNames, header))
                    .flatMap(Utils::stream)                                     // size == total # of contigs' from all successful assemblies
                    .filter(contig -> !contig.alignmentIntervals.isEmpty())     // filter out unmapped and contigs without primary alignments
                    .collect(Collectors.toList());
        }

        /**
         * Work on "successful" assembly and turn its contigs' alignments to custom {@link AlignedAssembly.AlignmentInterval} format.
         */
        @VisibleForTesting
        public static Iterable<AlignedAssembly.AlignedContig> forEachAssemblyNotExcuse(final AlignedAssemblyOrExcuse alignedAssembly,
                                                                                       final List<String> refNames,
                                                                                       final SAMFileHeader header) {

            final FermiLiteAssembly assembly = alignedAssembly.getAssembly();

            final List<List<BwaMemAlignment>> allAlignments = alignedAssembly.getContigAlignments();

            return IntStream.range(0, assembly.getNContigs())
                    .mapToObj( contigIdx -> {
                        final byte[] contigSequence = assembly.getContig(contigIdx).getSequence();
                        final String contigName = AlignedAssemblyOrExcuse.formatContigName(alignedAssembly.getAssemblyId(), contigIdx);
                        final List<AlignedAssembly.AlignmentInterval> arOfAContig = forEachContig(contigName, contigSequence, allAlignments.get(contigIdx), refNames, header);
                        return new AlignedAssembly.AlignedContig(contigName, contigSequence, arOfAContig);
                    } ).collect(Collectors.toList());
        }

        /**
         * Converts alignment records of the contig pointed to by {@code contigIdx} in a {@link FermiLiteAssembly} to custom {@link AlignedAssembly.AlignmentInterval} format.
         */
        @VisibleForTesting
        private static List<AlignedAssembly.AlignmentInterval> forEachContig(final String contigName, final byte[] contigSequence, final List<BwaMemAlignment> contigAlignments,
                                                                             final List<String> refNames, final SAMFileHeader header) {

            return contigAlignments.stream()
                    .filter( bwaMemAlignment ->  bwaMemAlignment.getRefId() >= 0 && SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwaMemAlignment.getSamFlag())) // mapped and not XA (i.e. not secondary)
                    .map(bwaMemAlignment -> BwaMemAlignmentUtils.applyAlignment(contigName, contigSequence, null, null, bwaMemAlignment, refNames, header, false, false))
                    .map(AlignedAssembly.AlignmentInterval::new)
                    .map(ar -> AssemblyAlignmentParser.breakGappedAlignment(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, contigSequence.length))
                    .flatMap(Utils::stream).collect(Collectors.toList());
        }
    }

}
