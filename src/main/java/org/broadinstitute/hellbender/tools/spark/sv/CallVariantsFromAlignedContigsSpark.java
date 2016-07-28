package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.collections4.IterableUtils;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.BreakpointAllele;
import org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.AlignAssembledContigsSpark.getContigsCollectionKeyedByBreakpointId;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = SparkProgramGroup.class)
public class CallVariantsFromAlignedContigsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "URL of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "Input file of contig alignments", shortName = "inputAlignments",
            fullName = "inputAlignments", optional = false)
    private String inputAlignments;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputAssemblies",
            fullName = "inputAssemblies", optional = false)
    private String inputAssemblies;

    // This class requires a reference parameter in 2bit format (to broadcast) and a reference in FASTA format
    // (to get a good sequence dictionary).
    // todo: document this better
    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        Broadcast<ReferenceMultiSource> broadcastReference = ctx.broadcast(getReference());

        final JavaRDD<AlignmentRegion> inputAlignedBreakpoints = ctx.textFile(inputAlignments).map(AlignAssembledContigsSpark::parseAlignedAssembledContigLine);

        final JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> alignmentRegionsKeyedByBreakpointAndContig = inputAlignedBreakpoints.mapToPair(alignmentRegion -> new Tuple2<>(new Tuple2<>(alignmentRegion.breakpointId, alignmentRegion.contigId), alignmentRegion)).groupByKey();

        final JavaPairRDD<String, ContigsCollection> breakpointIdsToContigsCollection = getContigsCollectionKeyedByBreakpointId(ctx, inputAssemblies).cache();

        final JavaPairRDD<Tuple2<String, String>, byte[]> contigSequences = breakpointIdsToContigsCollection.flatMapToPair(breakpointIdAndContigsCollection -> {
            final String breakpointId = breakpointIdAndContigsCollection._1;
            final ContigsCollection contigsCollection = breakpointIdAndContigsCollection._2;
            final List<Tuple2<Tuple2<String, String>, byte[]>> contigSequencesKeyedByBreakpoiuntIdAndContigId = contigsCollection.getContents().stream().map(pair -> new Tuple2<>(new Tuple2<>(breakpointId, pair._1.toString()), pair._2.toString().getBytes())).collect(Collectors.toList());
            return contigSequencesKeyedByBreakpoiuntIdAndContigId;
        });

        final JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> alignmentRegionsWithContigSequences = alignmentRegionsKeyedByBreakpointAndContig.join(contigSequences);

        JavaPairRDD<Tuple2<String, String>, AssembledBreakpoint> assembledBreakpointsByBreakpointIdAndContigId = alignmentRegionsWithContigSequences.flatMapValues(alignmentRegionsAndSequences -> CallVariantsFromAlignedContigsSpark.assembledBreakpointsFromAlignmentRegions(alignmentRegionsAndSequences._2, alignmentRegionsAndSequences._1));

        final JavaPairRDD<BreakpointAllele, Tuple2<Tuple2<String,String>, AssembledBreakpoint>> assembled3To5BreakpointsKeyedByPosition =
                assembledBreakpointsByBreakpointIdAndContigId
                        .mapToPair(CallVariantsFromAlignedContigsSpark::keyByBreakpointAllele)
                        .filter(CallVariantsFromAlignedContigsSpark::inversionBreakpointFilter);

        final JavaPairRDD<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>>> groupedBreakpoints = assembled3To5BreakpointsKeyedByPosition.groupByKey();

        final JavaRDD<VariantContext> variantContexts = groupedBreakpoints.map(breakpoints -> filterBreakpointsAndProduceVariants(breakpoints, broadcastReference)).cache();

        writeVariants(fastaReference, logger, variantContexts, getAuthenticatedGCSOptions(), outputPath);

    }

    public static void writeVariants(final String fastaReference, final Logger logger, final JavaRDD<VariantContext> variantContexts, final PipelineOptions pipelineOptions, final String outputPath) {

        final List<VariantContext> variants = variantContexts.collect();
        final List<VariantContext> sortedVariantsList = new ArrayList<>(variants);

        final ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(pipelineOptions, fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final SAMSequenceDictionary referenceSequenceDictionary = referenceMultiSource.getReferenceSequenceDictionary(null);

        sortedVariantsList.sort((VariantContext v1, VariantContext v2) -> IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary));

        logger.info("Called " + variants.size() + " candidate inversions");
        final VCFHeader header = getVcfHeader(referenceSequenceDictionary);

        writeVariants(outputPath, sortedVariantsList, pipelineOptions, "inversions.vcf", header, referenceSequenceDictionary);

        final List<VariantContext> hqmappingVariants = sortedVariantsList.stream().filter(v -> v.getAttributeAsInt(GATKSVVCFHeaderLines.HQ_MAPPINGS, 0) > 0).collect(Collectors.toList());
        logger.info("Called " + hqmappingVariants.size() + " high-quality inversions");
        writeVariants(outputPath, hqmappingVariants, pipelineOptions, "hq_inversions.vcf", header, referenceSequenceDictionary);
    }

    private static VCFHeader getVcfHeader(final SAMSequenceDictionary referenceSequenceDictionary) {
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(referenceSequenceDictionary);
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        GATKSVVCFHeaderLines.vcfHeaderLines.values().forEach(header::addMetaDataLine);
        return header;
    }

    private static void writeVariants(final String outputPath, final List<VariantContext> variantsArrayList, final PipelineOptions pipelineOptions, final String fileName, final VCFHeader header, final SAMSequenceDictionary referenceSequenceDictionary) {
        try (final OutputStream outputStream = new BufferedOutputStream(
                BucketUtils.createFile(outputPath + "/" + fileName, pipelineOptions))) {

            final VariantContextWriter vcfWriter = getVariantContextWriter(outputStream, referenceSequenceDictionary);

            vcfWriter.writeHeader(header);
            variantsArrayList.forEach(vcfWriter::add);
            vcfWriter.close();

        } catch (IOException e) {
            throw new GATKException("Could not create output file", e);
        }
    }

    private static VariantContextWriter getVariantContextWriter(final OutputStream outputStream, final SAMSequenceDictionary referenceSequenceDictionary) {
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                                                            .clearOptions()
                                                            .setOutputStream(outputStream);

        if (null != referenceSequenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceSequenceDictionary);
        }
        // todo: remove this when things are solid?
        vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        for (Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

    public static Iterable<AssembledBreakpoint> assembledBreakpointsFromAlignmentRegions(final byte[] contigSequence, final Iterable<AlignmentRegion> alignmentRegionsIterable) {
        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionsIterable);
        if (alignmentRegions.size() > 1) {
            alignmentRegions.sort(Comparator.comparing(a -> a.startInAssembledContig));
        }
        return ContigAligner.getAssembledBreakpointsFromAlignmentRegions(contigSequence, alignmentRegions);
    }

    @VisibleForTesting
    public static VariantContext filterBreakpointsAndProduceVariants(final Tuple2<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>>> assembledBreakpointsPerAllele, final Broadcast<ReferenceMultiSource> broadcastReference) throws IOException {
        int numAssembledBreakpoints = 0;
        int highMqMappings = 0;
        int midMqMappings = 0;
        int lowMqMappings = 0;
        int maxAlignLength = 0;

        final Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>> assembledBreakpoints = assembledBreakpointsPerAllele._2;
        final List<Tuple2<Tuple2<String, String>, AssembledBreakpoint>> assembledBreakpointsList = IterableUtils.toList(assembledBreakpoints);
        final int numBreakpoints = assembledBreakpointsList.size();
        final List<Integer> mqs = new ArrayList<>(numBreakpoints);
        final List<Integer> alignLengths = new ArrayList<>(numBreakpoints);
        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final List<String> breakpointIds = new ArrayList<>(numBreakpoints);
        final List<String> assembledContigIds = new ArrayList<>(numBreakpoints);

        for (final Tuple2<Tuple2<String,String>, AssembledBreakpoint> assembledBreakpointPair : assembledBreakpointsList) {
            final AssembledBreakpoint assembledBreakpoint = assembledBreakpointPair._2;
            numAssembledBreakpoints = numAssembledBreakpoints + 1;
            final int assembledBreakpointMapq = Math.min(assembledBreakpoint.region1.mqual, assembledBreakpoint.region2.mqual);
            if (assembledBreakpointMapq == 60) {
                highMqMappings = highMqMappings + 1;
            } else if (assembledBreakpointMapq > 0) {
                midMqMappings = midMqMappings + 1;
            } else {
                lowMqMappings = lowMqMappings + 1;
            }
            mqs.add(assembledBreakpointMapq);
            final int assembledBreakpointAlignmentLength = Math.min(assembledBreakpoint.region1.referenceInterval.size(), assembledBreakpoint.region2.referenceInterval.size());
            alignLengths.add(assembledBreakpointAlignmentLength);
            maxAlignLength = Math.max(maxAlignLength, assembledBreakpointAlignmentLength);
            breakpointIds.add(assembledBreakpointPair._1._1);
            assembledContigIds.add(assembledBreakpoint.contigId);
        }

        return createVariant(numAssembledBreakpoints, highMqMappings, mqs, alignLengths, maxAlignLength, breakpointAllele, breakpointIds, assembledContigIds, broadcastReference.getValue());
    }

    @VisibleForTesting
    private static VariantContext createVariant(final int numAssembledBreakpoints, final int highMqMappings, final List<Integer> mqs, final List<Integer> alignLengths, final int maxAlignLength, final BreakpointAllele breakpointAllele, final List<String> breakpointIds, final List<String> assembledContigIds, final ReferenceMultiSource reference) throws IOException {
        final String contig = breakpointAllele.leftAlignedLeftBreakpoint.getContig();
        final int start = breakpointAllele.leftAlignedLeftBreakpoint.getStart();
        final int end = breakpointAllele.leftAlignedRightBreakpoint.getStart();

        final Allele refAllele = Allele.create(new String(reference.getReferenceBases(null, new SimpleInterval(contig, start, start)).getBases()), true);
        final Allele altAllele = Allele.create("<INV>");
        final List<Allele> vcAlleles = new ArrayList<>(2);
        vcAlleles.add(refAllele);
        vcAlleles.add(altAllele);

        VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(contig)
                .start(start)
                .stop(end)
                .id(getInversionId(breakpointAllele))
                .alleles(vcAlleles)
                .attribute(VCFConstants.END_KEY, end)
                .attribute(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, end - start)
                .attribute(GATKSVVCFHeaderLines.TOTAL_MAPPINGS, numAssembledBreakpoints)
                .attribute(GATKSVVCFHeaderLines.HQ_MAPPINGS, highMqMappings)
                .attribute(GATKSVVCFHeaderLines.MAPPING_QUALITIES, mqs.stream().map(String::valueOf).collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.ALIGN_LENGTHS, alignLengths.stream().map(String::valueOf).collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, maxAlignLength)
                .attribute(GATKSVVCFHeaderLines.BREAKPOINT_IDS, breakpointIds.stream().collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, assembledContigIds.stream().map(s -> s.replace(" ", "_")).collect(Collectors.joining(",")));

        if (breakpointAllele.insertedSequence.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, breakpointAllele.insertedSequence);
        }

        if (is5To3Inversion(breakpointAllele)) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_5_TO_3, "");
        }

        if (is3To5Inversion(breakpointAllele)) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }


        return vcBuilder.make();
    }

    private static boolean is3To5Inversion(final BreakpointAllele breakpointAllele) {
        return ! breakpointAllele.fiveToThree && breakpointAllele.threeToFive;
    }

    private static boolean is5To3Inversion(final BreakpointAllele breakpointAllele) {
        return breakpointAllele.fiveToThree && ! breakpointAllele.threeToFive;
    }

    private static String getInversionId(final BreakpointAllele breakpointAllele) {
        final String invType;
        if (is5To3Inversion(breakpointAllele)) {
            invType = "INV_5_TO_3";
        } else if (is3To5Inversion(breakpointAllele)) {
            invType = "INV_3_TO_5";
        } else {
            // I don't think this should happen
            invType = "INV_UNKNOWN";
        }

        return invType + "_" + breakpointAllele.leftAlignedLeftBreakpoint.getContig() + "_" + breakpointAllele.leftAlignedLeftBreakpoint.getStart() + "_" + breakpointAllele.leftAlignedRightBreakpoint.getStart();
    }


    public static Tuple2<BreakpointAllele, Tuple2<Tuple2<String,String>, AssembledBreakpoint>> keyByBreakpointAllele(final Tuple2<Tuple2<String,String>, AssembledBreakpoint> breakpointIdAndAssembledBreakpoint) {
        final BreakpointAllele breakpointAllele = breakpointIdAndAssembledBreakpoint._2.getBreakpointAllele();
        return new Tuple2<>(breakpointAllele, new Tuple2<>(breakpointIdAndAssembledBreakpoint._1, breakpointIdAndAssembledBreakpoint._2));
    }

    public static Boolean inversionBreakpointFilter(final Tuple2<BreakpointAllele, Tuple2<Tuple2<String,String>,AssembledBreakpoint>> breakpointAlleleTuple2Tuple2) {
        final BreakpointAllele breakpointAllele = breakpointAlleleTuple2Tuple2._1;
        return breakpointAllele.leftAlignedLeftBreakpoint.getContig().equals(breakpointAllele.leftAlignedRightBreakpoint.getContig()) && (breakpointAllele.fiveToThree || breakpointAllele.threeToFive);
    }
}
