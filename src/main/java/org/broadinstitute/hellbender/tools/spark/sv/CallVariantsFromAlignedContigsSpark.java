package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.IterableUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.BreakpointAllele;
import org.broadinstitute.hellbender.tools.spark.sv.RunSGAViaProcessBuilderOnSpark.ContigsCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
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

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

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
                        .filter(CallVariantsFromAlignedContigsSpark::inversionBreakpointFilter)
                        .mapToPair(CallVariantsFromAlignedContigsSpark::keyByBreakpointAllele);

        final JavaPairRDD<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>>> groupedBreakpoints = assembled3To5BreakpointsKeyedByPosition.groupByKey();

        final JavaRDD<VariantContext> variantContexts = groupedBreakpoints.map(CallVariantsFromAlignedContigsSpark::filterBreakpointsAndProduceVariants).cache();

        final List<VariantContext> variants = variantContexts.collect();
        final List<VariantContext> variantsArrayList = new ArrayList<>(variants);
        variantsArrayList.sort((VariantContext v1, VariantContext v2) -> IntervalUtils.compareLocatables(v1, v2, getReferenceSequenceDictionary()));

        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();

        final VCFHeader header = getVcfHeader();

        writeVariants(variantsArrayList, pipelineOptions, "inversions.vcf", header);

        final List<VariantContext> hqmappingVariants = variantsArrayList.stream().filter(v -> v.getAttributeAsInt(GATKSVVCFHeaderLines.HQ_MAPPINGS, 0) > 0).collect(Collectors.toList());
        writeVariants(hqmappingVariants, pipelineOptions, "hq_inversions.vcf", header);

    }

    private VCFHeader getVcfHeader() {
        final VCFHeader header = new VCFHeader();
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        GATKSVVCFHeaderLines.vcfHeaderLines.values().forEach(header::addMetaDataLine);
        return header;
    }

    private void writeVariants(final List<VariantContext> variantsArrayList, final PipelineOptions pipelineOptions, final String fileName, final VCFHeader header) {
        try (final OutputStream outputStream = new BufferedOutputStream(
                BucketUtils.createFile(outputPath  + "/" + fileName, pipelineOptions))) {

            final VariantContextWriter vcfWriter = getVariantContextWriter(outputStream);

            vcfWriter.writeHeader(header);
            variantsArrayList.forEach(vcfWriter::add);
            vcfWriter.close();

        } catch (IOException e) {
            throw new GATKException("Could not create output file", e);
        }
    }

    private VariantContextWriter getVariantContextWriter(final OutputStream outputStream) {
        final SAMSequenceDictionary referenceDictionary = getReferenceSequenceDictionary();
        VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder()
                                                            .clearOptions()
                                                            .setOutputStream(outputStream)
                                                            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);

        if (null != referenceDictionary) {
            vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
        }
        vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        for (Options opt : new Options[]{}) {
            vcWriterBuilder = vcWriterBuilder.setOption(opt);
        }

        return vcWriterBuilder.build();
    }

    private static Iterable<AssembledBreakpoint> assembledBreakpointsFromAlignmentRegions(final byte[] contigSequence, final Iterable<AlignmentRegion> alignmentRegionsIterable) {
        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionsIterable);
        alignmentRegions.sort(Comparator.comparing(a -> a.startInAssembledContig));
        return ContigAligner.getAssembledBreakpointsFromAlignmentRegions(contigSequence, alignmentRegions);
    }

    @VisibleForTesting
    private static VariantContext filterBreakpointsAndProduceVariants(final Tuple2<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>>> assembledBreakpointsPerAllele) throws IOException {
        int numAssembledBreakpoints = 0;
        int highMqMappings = 0;
        int midMqMappings = 0;
        int lowMqMappings = 0;
        int maxAlignLength = 0;
        final List<Integer> mqs = new ArrayList<>(10);
        final List<Integer> alignLengths = new ArrayList<>(10);
        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final List<String> breakpointIds = new ArrayList<>(10);
        final List<String> assembledContigIds = new ArrayList<>(10);
        final Iterable<Tuple2<Tuple2<String,String>, AssembledBreakpoint>> assembledBreakpoints = assembledBreakpointsPerAllele._2;
        for (final Tuple2<Tuple2<String,String>, AssembledBreakpoint> assembledBreakpointPair : assembledBreakpoints) {
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

        return createVariant(numAssembledBreakpoints, highMqMappings, mqs, alignLengths, maxAlignLength, breakpointAllele, breakpointIds, assembledContigIds);
    }

    @VisibleForTesting
    private static VariantContext createVariant(final int numAssembledBreakpoints, final int highMqMappings, final List<Integer> mqs, final List<Integer> alignLengths, final int maxAlignLength, final BreakpointAllele breakpointAllele, final List<String> breakpointIds, final List<String> assembledContigIds) throws IOException {
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("<INV>");
        final List<Allele> vcAlleles = new ArrayList<>(2);
        vcAlleles.add(refAllele);
        vcAlleles.add(altAllele);

        VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(breakpointAllele.leftAlignedLeftBreakpoint.getContig())
                .start(breakpointAllele.leftAlignedLeftBreakpoint.getStart())
                .stop(breakpointAllele.leftAlignedRightBreakpoint.getStart())
                .id(getInversionId(breakpointAllele))
                .alleles(vcAlleles)
                .attribute(VCFConstants.END_KEY, breakpointAllele.leftAlignedRightBreakpoint.getStart())
                .attribute(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.toString())
                .attribute(GATKSVVCFHeaderLines.SVLEN, breakpointAllele.leftAlignedRightBreakpoint.getStart() - breakpointAllele.leftAlignedLeftBreakpoint.getStart())
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
        return ! breakpointAllele.left5Prime && breakpointAllele.right5Prime;
    }

    private static boolean is5To3Inversion(final BreakpointAllele breakpointAllele) {
        return breakpointAllele.left5Prime && ! breakpointAllele.right5Prime;
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

    @VisibleForTesting
    static boolean inversionBreakpointFilter(final Tuple2<Tuple2<String, String>, AssembledBreakpoint> assembledBreakpoint) {
        final AlignmentRegion region1 = assembledBreakpoint._2.region1;
        final AlignmentRegion region2 = assembledBreakpoint._2.region2;
        return region1.referenceInterval.getContig().equals(region2.referenceInterval.getContig()) &&
                (region1.forwardStrand && ! region2.forwardStrand || ! region1.forwardStrand && region2.forwardStrand) &&
                ! region1.referenceInterval.overlaps(region2.referenceInterval);
    }

    private static Tuple2<BreakpointAllele, Tuple2<Tuple2<String,String>, AssembledBreakpoint>> keyByBreakpointAllele(final Tuple2<Tuple2<String,String>, AssembledBreakpoint> breakpointIdAndAssembledBreakpoint) {
        final BreakpointAllele breakpointAllele = breakpointIdAndAssembledBreakpoint._2.getBreakpointAllele();
        return new Tuple2<>(breakpointAllele, new Tuple2<>(breakpointIdAndAssembledBreakpoint._1, breakpointIdAndAssembledBreakpoint._2));
    }
}
