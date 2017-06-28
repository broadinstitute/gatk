package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import scala.Serializable;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tool to run the sv pipeline up to and including variant discovery
 */
@CommandLineProgramProperties(summary="Master tool to run the structural variation discovery pipeline",
        oneLineSummary="Master tool to run the structural variation discovery pipeline",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
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
    private StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection evidenceAndAssemblyArgs
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
                        .gatherEvidenceAndWriteContigSamFile(ctx, evidenceAndAssemblyArgs, header, getUnfilteredReads(),
                                outputSAM, localLogger);
        if (alignedAssemblyOrExcuseList.isEmpty()) return;

        // parse the contig alignments and extract necessary information
        @SuppressWarnings("unchecked")
        final JavaRDD<AlignedContig> parsedAlignments = new InMemoryAlignmentParser(ctx, alignedAssemblyOrExcuseList, header, localLogger).getAlignedContigs();
        if(parsedAlignments.isEmpty()) return;

        // discover variants and write to vcf
        DiscoverVariantsFromContigAlignmentsSAMSpark
                .discoverVariantsAndWriteVCF(parsedAlignments, discoverStageArgs.fastaReference,
                        ctx.broadcast(getReference()), pipelineOptions, vcfOutputFileName, localLogger);
    }

    public static final class InMemoryAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaSparkContext ctx;
        private final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        private final SAMFileHeader header;
        private final Logger toolLogger;


        InMemoryAlignmentParser(final JavaSparkContext ctx,
                                final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                final SAMFileHeader header,
                                final Logger toolLogger) {
            this.ctx = ctx;
            this.alignedAssemblyOrExcuseList = alignedAssemblyOrExcuseList;
            this.header = header;
            this.toolLogger = toolLogger;
        }

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {

            // here we have two options, one is going through the route "BwaMemAlignment -> SAM -> GATKRead -> SAM -> AlignmentInterval"
            //                           which is the route if the discovery pipeline is run by "FindBreakpointEvidenceSpark -> write sam file -> load sam file -> DiscoverVariantsFromContigAlignmentsSAMSpark"
            //                         , the other is to go directly "BwaMemAlignment -> AlignmentInterval" by calling into {@code filterAndConvertToAlignedContigDirect()}, which is faster but not used here.
            //                         ; the two routes are tested to be generating the same output via {@code AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
            return filterAndConvertToAlignedContigViaSAM(alignedAssemblyOrExcuseList, header, ctx, toolLogger);
        }

        /**
         * Filters out "failed" assemblies, unmapped and secondary (i.e. "XA") alignments, and
         * turn the alignments of contigs into custom {@link AlignedAssembly.AlignmentInterval} format.
         * Should be generating the same output as {@link #filterAndConvertToAlignedContigDirect(Iterable, List, SAMFileHeader)};
         * and currently this assertion is tested in {@see AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
         */
        @VisibleForTesting
        public static JavaRDD<AlignedContig> filterAndConvertToAlignedContigViaSAM(final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                                                                   final SAMFileHeader header,
                                                                                   final JavaSparkContext ctx,
                                                                                   final Logger toolLogger) {

            final SAMFileHeader cleanHeader = new SAMFileHeader(header.getSequenceDictionary());
            final List<String> refNames = AlignedAssemblyOrExcuse.getRefNames(header);

            return ctx.parallelize(alignedAssemblyOrExcuseList)
                    .filter(AlignedAssemblyOrExcuse::isNotFailure)
                    .flatMap(alignedAssemblyNoExcuse -> {
                                final FermiLiteAssembly assembly = alignedAssemblyNoExcuse.getAssembly();
                                final int assemblyId = alignedAssemblyNoExcuse.getAssemblyId();
                                final List<List<BwaMemAlignment>> allAlignmentsOfThisAssembly = alignedAssemblyNoExcuse.getContigAlignments();
                                final int nContigs = assembly.getNContigs();
                                return IntStream.range(0, nContigs)
                                        .mapToObj(contigIdx ->
                                                AlignedAssemblyOrExcuse.toSAMStreamForOneContig(cleanHeader, refNames, assemblyId, contigIdx,
                                                        assembly.getContig(contigIdx).getSequence(), allAlignmentsOfThisAssembly.get(contigIdx))
                                        ).iterator();
                            }
                    )
                    .map(forOneContig -> forOneContig.filter(sam -> !sam.getReadUnmappedFlag() && !sam.getNotPrimaryAlignmentFlag()).collect(Collectors.toList()))
                    .filter(list -> !list.isEmpty()) // not filtering on the stream directly to avoid consuming the stream so next operations would not throw
                    .map(forOneContig ->
                            DiscoverVariantsFromContigAlignmentsSAMSpark.
                                    SAMFormattedContigAlignmentParser.
                                    parseReadsAndBreakGaps(forOneContig,
                                            header, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, toolLogger ));
        }

        /**
         * Filter out "failed" assemblies, unmapped and secondary (i.e. "XA") alignments, and
         * turn the alignments of contigs into custom {@link AlignedAssembly.AlignmentInterval} format.
         * Should be generating the same output as {@link #filterAndConvertToAlignedContigViaSAM(List, SAMFileHeader, JavaSparkContext, Logger)};
         * and currently this assertion is tested in {@see AlignedContigGeneratorUnitTest#testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute()}
         */
        @VisibleForTesting
        public static List<AlignedContig> filterAndConvertToAlignedContigDirect(final Iterable<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseIterable,
                                                                                final List<String> refNames, final SAMFileHeader header) {

            return Utils.stream(alignedAssemblyOrExcuseIterable)
                    .filter(AlignedAssemblyOrExcuse::isNotFailure)
                    .map(alignedAssembly -> getAlignedContigsInOneAssembly(alignedAssembly, refNames, header))
                    .flatMap(Utils::stream)                                     // size == total # of contigs' from all successful assemblies
                    .filter(contig -> !contig.alignmentIntervals.isEmpty())     // filter out unmapped and contigs without primary alignments
                    .collect(Collectors.toList());
        }

        /**
         * Work on "successful" assembly and turn its contigs' alignments to custom {@link AlignedAssembly.AlignmentInterval} format.
         */
        @VisibleForTesting
        public static Iterable<AlignedContig> getAlignedContigsInOneAssembly(final AlignedAssemblyOrExcuse alignedAssembly,
                                                                             final List<String> refNames,
                                                                             final SAMFileHeader header) {

            final FermiLiteAssembly assembly = alignedAssembly.getAssembly();

            final List<List<BwaMemAlignment>> allAlignments = alignedAssembly.getContigAlignments();

            return IntStream.range(0, assembly.getNContigs())
                    .mapToObj( contigIdx -> {
                        final byte[] contigSequence = assembly.getContig(contigIdx).getSequence();
                        final String contigName = AlignedAssemblyOrExcuse.formatContigName(alignedAssembly.getAssemblyId(), contigIdx);
                        final List<AlignedAssembly.AlignmentInterval> arOfAContig = getAlignmentsForOneContig(contigName, contigSequence, allAlignments.get(contigIdx), refNames, header);
                        return new AlignedContig(contigName, contigSequence, arOfAContig);
                    } ).collect(Collectors.toList());
        }

        /**
         * Converts alignment records of the contig pointed to by {@code contigIdx} in a {@link FermiLiteAssembly} to custom {@link AlignedAssembly.AlignmentInterval} format.
         * Filters out unmapped and ambiguous mappings (i.e. XA).
         */
        @VisibleForTesting
        private static List<AlignedAssembly.AlignmentInterval> getAlignmentsForOneContig(final String contigName,
                                                                                         final byte[] contigSequence,
                                                                                         final List<BwaMemAlignment> contigAlignments,
                                                                                         final List<String> refNames,
                                                                                         final SAMFileHeader header) {

            return contigAlignments.stream()
                    .filter( bwaMemAlignment ->  bwaMemAlignment.getRefId() >= 0 && SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(bwaMemAlignment.getSamFlag())) // mapped and not XA (i.e. not secondary)
                    .map(bwaMemAlignment -> BwaMemAlignmentUtils.applyAlignment(contigName, contigSequence, null,
                            null, bwaMemAlignment, refNames, header, false, false))
                    .map(AlignedAssembly.AlignmentInterval::new)
                    .map(ar -> GappedAlignmentSplitter.split(ar, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, contigSequence.length))
                    .flatMap(Utils::stream).collect(Collectors.toList());
        }
    }

}
