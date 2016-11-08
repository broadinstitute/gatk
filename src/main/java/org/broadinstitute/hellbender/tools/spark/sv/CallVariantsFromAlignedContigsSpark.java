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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.ContigsCollection.loadContigsCollectionKeyedByAssemblyId;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class CallVariantsFromAlignedContigsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    public static final Integer DEFAULT_MIN_ALIGNMENT_LENGTH = 50;

    private static final Logger log = LogManager.getLogger(CallVariantsFromAlignedContigsSpark.class);

    public static final String INVERSIONS_OUTPUT_VCF = "inversions.vcf";

    @Argument(doc = "URI of the output path", shortName = "outputPath",
            fullName = "outputPath", optional = false)
    private String outputPath;

    @Argument(doc = "Input file of contig alignments", shortName = "inputAlignments",
            fullName = "inputAlignments", optional = false)
    private String inputAlignments;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputAssemblies",
            fullName = "inputAssemblies", optional = false)
    private String inputAssemblies;

    @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
            fullName = "minAlignLength", optional = true)
    private Integer minAlignLength = CallVariantsFromAlignedContigsSpark.DEFAULT_MIN_ALIGNMENT_LENGTH;

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

        final JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> alignmentRegionsWithContigSequences = prepAlignmentRegionsForCalling(ctx);

        final Integer minAlignLengthFinal = this.minAlignLength;
        callVariantsFromAlignmentRegionsAndWriteVariants(broadcastReference, alignmentRegionsWithContigSequences, minAlignLengthFinal, fastaReference, getAuthenticatedGCSOptions(), outputPath);


    }

    /**
     * Loads the alignment regions from the text file they are in; converts them to a PairRDD keyed by breakpoint and contig ID;
     * loads the assembled contigs for all assemblies and uses them to add the contig sequence to each item in the PairRDD.
     */
    private JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> prepAlignmentRegionsForCalling(final JavaSparkContext ctx) {
        final JavaRDD<AlignmentRegion> inputAlignedContigs = ctx.textFile(inputAlignments).map(ContigsCollection::parseAlignedAssembledContigLine);

        final JavaPairRDD<Tuple2<String, String>, Iterable<AlignmentRegion>> alignmentRegionsKeyedByBreakpointAndContig = inputAlignedContigs.mapToPair(alignmentRegion -> new Tuple2<>(new Tuple2<>(alignmentRegion.assemblyId, alignmentRegion.contigId), alignmentRegion)).groupByKey();

        final JavaPairRDD<String, ContigsCollection> assemblyIdsToContigCollections = loadContigsCollectionKeyedByAssemblyId(ctx, inputAssemblies);

        final JavaPairRDD<Tuple2<String, String>, byte[]> contigSequences = assemblyIdsToContigCollections.flatMapToPair(assemblyIdAndContigsCollection -> {
            final String assemblyId = assemblyIdAndContigsCollection._1;
            final ContigsCollection contigsCollection = assemblyIdAndContigsCollection._2;
            return contigsCollection.getContents().stream().map(pair -> new Tuple2<>(new Tuple2<>(assemblyId, pair._1.toString()), pair._2.toString().getBytes())).collect(Collectors.toList()).iterator();
        });

        return alignmentRegionsKeyedByBreakpointAndContig.join(contigSequences);
    }

    /**
     * This method processes an RDD containing alignment regions, scanning for split alignments which match a set of filtering
     * criteria, and emitting a list of VariantContexts representing SVs for split alignments that pass. It then writes the variants
     * to the output path.
     *
     * The input RDD is of the form:
     *
     * Key: Tuple2 of two Strings: the assembly ID and the contig ID that the alignments come from
     * Value: Tuple2 of:
     *     {@code Iterable<AlignmentRegion>} AlignmentRegion objects representing all alignments for the contig
     *     A byte array with the sequence content of the contig
     *
     * FASTA and Broadcast references are both required because 2bit Broadcast references currently order their
     * sequence dictionaries in a scrambled order, see https://github.com/broadinstitute/gatk/issues/2037.
     *
     * @param broadcastReference The broadcast handle to the reference (used to populate reference bases)
     * @param alignmentRegionsWithContigSequences A data structure as described above, where a list of AlignmentRegions and the sequence of the contig are keyed by a tuple of Assembly ID and Contig ID
     * @param minAlignLength The minimum length of alignment regions flanking the breakpoint to emit an SV variant
     * @param fastaReference The reference in FASTA format, used to get a sequence dictionary and sort the variants according to it
     * @param pipelineOptions GCS pipeline option for creating and writing output files
     * @param outputPath Path to write the output files
     * @return An RDD of VariantContexts representing SVs called from breakpoint alignments
     */
    protected static void callVariantsFromAlignmentRegionsAndWriteVariants(final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                           final JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> alignmentRegionsWithContigSequences,
                                                                           final Integer minAlignLength,
                                                                           final String fastaReference,
                                                                           final PipelineOptions pipelineOptions,
                                                                           final String outputPath) {
        final JavaRDD<VariantContext> variants = callVariantsFromAlignmentRegions(broadcastReference, alignmentRegionsWithContigSequences, minAlignLength);
        writeVariants(fastaReference, variants, pipelineOptions, outputPath);
    }

    /**
     * This method processes an RDD containing alignment regions, scanning for split alignments which match a set of filtering
     * criteria, and emitting a list of VariantContexts representing SVs for split alignments that pass.
     *
     * The input RDD is of the form:
     *
     * Key: Tuple2 of two Strings: the assembly ID and the contig ID that the alignments come from
     * Value: Tuple2 of:
     *     {@code Iterable<AlignmentRegion>} AlignmentRegion objects representing all alignments for the contig
     *     A byte array with the sequence content of the contig
     *
     * @param broadcastReference The broadcast handle to the reference (used to populate reference bases)
     * @param alignmentRegionsWithContigSequences A data structure as described above, where a list of AlignmentRegions and the sequence of the contig are keyed by a tuple of Assembly ID and Contig ID
     * @param minAlignLength The minimum length of alignment regions flanking the breakpoint to emit an SV variant
     * @return An RDD of VariantContexts representing SVs called from breakpoint alignments
     */
    protected static JavaRDD<VariantContext> callVariantsFromAlignmentRegions(final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                              final JavaPairRDD<Tuple2<String, String>, Tuple2<Iterable<AlignmentRegion>, byte[]>> alignmentRegionsWithContigSequences,
                                                                              final Integer minAlignLength) {
        JavaPairRDD<Tuple2<String, String>, BreakpointAlignment> assembledBreakpointsByBreakpointIdAndContigId =
                alignmentRegionsWithContigSequences.flatMapValues(alignmentRegionsAndSequences ->
                        CallVariantsFromAlignedContigsSpark.assembledBreakpointsFromAlignmentRegions(alignmentRegionsAndSequences._2, alignmentRegionsAndSequences._1, minAlignLength));

        final JavaPairRDD<BreakpointAllele, Tuple2<Tuple2<String,String>, BreakpointAlignment>> inversionBreakpointAlignmentsKeyedByAllele =
                assembledBreakpointsByBreakpointIdAndContigId
                        .mapToPair(CallVariantsFromAlignedContigsSpark::keyByBreakpointAllele)
                        .filter(pair -> pair._1().isInversion());

        final JavaPairRDD<BreakpointAllele, Iterable<Tuple2<Tuple2<String,String>, BreakpointAlignment>>> groupedBreakpoints =
                inversionBreakpointAlignmentsKeyedByAllele.groupByKey();

        return groupedBreakpoints.map(breakpoints -> getVariantContextForBreakpointAlleleAlignmentList(breakpoints, broadcastReference));
    }

    static List<BreakpointAlignment> getBreakpointAlignmentsFromAlignmentRegions(final byte[] sequence, final List<AlignmentRegion> alignmentRegionList, final Integer minAlignLength) {
        if (alignmentRegionList.isEmpty()) {
            return new ArrayList<>();
        }
        final List<BreakpointAlignment> results = new ArrayList<>(alignmentRegionList.size() - 1);
        final Iterator<AlignmentRegion> iterator = alignmentRegionList.iterator();
        final List<String> insertionAlignmentRegions = new ArrayList<>();
        if ( iterator.hasNext() ) {
            AlignmentRegion current = iterator.next();
            while (treatAlignmentRegionAsInsertion(current) && iterator.hasNext()) {
                current = iterator.next();
            }
            while ( iterator.hasNext() ) {
                final AlignmentRegion next = iterator.next();
                if (currentAlignmentRegionIsTooSmall(current, next, minAlignLength)) {
                    continue;
                }

                if (treatNextAlignmentRegionInPairAsInsertion(current, next, minAlignLength)) {
                    if (iterator.hasNext()) {
                        insertionAlignmentRegions.add(next.toPackedString());
                        continue;
                    } else {
                        break;
                    }
                }

                final AlignmentRegion previous = current;
                current = next;

                final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);

                String homology = getHomology(current, previous, sequenceCopy);
                String insertedSequence = getInsertedSequence(current, previous, sequenceCopy);

                final BreakpointAlignment breakpointAlignment = new BreakpointAlignment(current.contigId, previous, current, insertedSequence, homology, insertionAlignmentRegions);

                results.add(breakpointAlignment);
            }
        }
        return results;
    }

    private static String getInsertedSequence(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String insertedSequence = "";
        if (previous.endInAssembledContig < current.startInAssembledContig - 1) {

            final int insertionStart;
            final int insertionEnd;

            insertionStart = previous.endInAssembledContig + 1;
            insertionEnd = current.startInAssembledContig - 1;

            final byte[] insertedSequenceBytes = Arrays.copyOfRange(sequenceCopy, insertionStart - 1, insertionEnd);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(insertedSequenceBytes, 0, insertedSequenceBytes.length);
            }
            insertedSequence = new String(insertedSequenceBytes);
        }
        return insertedSequence;
    }

    private static String getHomology(final AlignmentRegion current, final AlignmentRegion previous, final byte[] sequenceCopy) {
        String homology = "";
        if (previous.endInAssembledContig >= current.startInAssembledContig) {
            final byte[] homologyBytes = Arrays.copyOfRange(sequenceCopy, current.startInAssembledContig - 1, previous.endInAssembledContig);
            if (previous.referenceInterval.getStart() > current.referenceInterval.getStart()) {
                SequenceUtil.reverseComplement(homologyBytes, 0, homologyBytes.length);
            }
            homology = new String(homologyBytes);
        }
        return homology;
    }

    private static boolean currentAlignmentRegionIsTooSmall(final AlignmentRegion current, final AlignmentRegion next, final Integer minAlignLength) {
        return current.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength;
    }

    @VisibleForTesting
    static boolean treatNextAlignmentRegionInPairAsInsertion(AlignmentRegion current, AlignmentRegion next, final Integer minAlignLength) {
        return treatAlignmentRegionAsInsertion(next) ||
                (next.referenceInterval.size() - current.overlapOnContig(next) < minAlignLength) ||
                current.referenceInterval.contains(next.referenceInterval) ||
                next.referenceInterval.contains(current.referenceInterval);
    }

    private static boolean treatAlignmentRegionAsInsertion(final AlignmentRegion next) {
        return next.mapqual < 60;
    }

    private static void writeVariants(final String fastaReference, final JavaRDD<VariantContext> variantContexts, final PipelineOptions pipelineOptions, final String outputPath) {

        final List<VariantContext> variants = variantContexts.collect();
        final List<VariantContext> sortedVariantsList = new ArrayList<>(variants);

        final ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(pipelineOptions, fastaReference, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final SAMSequenceDictionary referenceSequenceDictionary = referenceMultiSource.getReferenceSequenceDictionary(null);

        sortedVariantsList.sort((VariantContext v1, VariantContext v2) -> IntervalUtils.compareLocatables(v1, v2, referenceSequenceDictionary));

        log.info("Called " + variants.size() + " inversions");
        final VCFHeader header = getVcfHeader(referenceSequenceDictionary);

        writeVariants(outputPath, sortedVariantsList, pipelineOptions, INVERSIONS_OUTPUT_VCF, header, referenceSequenceDictionary);
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

    private static Iterable<BreakpointAlignment> assembledBreakpointsFromAlignmentRegions(final byte[] contigSequence, final Iterable<AlignmentRegion> alignmentRegionsIterable, final Integer minAlignLength) {
        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionsIterable);
        if (alignmentRegions.size() > 1) {
            alignmentRegions.sort(Comparator.comparing(a -> a.startInAssembledContig));
        }
        return getBreakpointAlignmentsFromAlignmentRegions(contigSequence, alignmentRegions, minAlignLength);
    }

    @VisibleForTesting
    private static VariantContext getVariantContextForBreakpointAlleleAlignmentList(final Tuple2<BreakpointAllele, Iterable<Tuple2<Tuple2<String, String>, BreakpointAlignment>>> assembledBreakpointsPerAllele, final Broadcast<ReferenceMultiSource> broadcastReference) throws IOException {
        int numAssembledBreakpoints = 0;
        int highMqMappings = 0;
        int maxAlignLength = 0;

        final Iterable<Tuple2<Tuple2<String,String>, BreakpointAlignment>> assembledBreakpoints = assembledBreakpointsPerAllele._2;
        final List<Tuple2<Tuple2<String, String>, BreakpointAlignment>> assembledBreakpointsList = IterableUtils.toList(assembledBreakpoints);
        final int numBreakpoints = assembledBreakpointsList.size();
        final List<Integer> mqs = new ArrayList<>(numBreakpoints);
        final List<Integer> alignLengths = new ArrayList<>(numBreakpoints);
        final BreakpointAllele breakpointAllele = assembledBreakpointsPerAllele._1;
        final List<String> breakpointIds = new ArrayList<>(numBreakpoints);
        final List<String> assembledContigIds = new ArrayList<>(numBreakpoints);

        final List<String> insertionMappings = new ArrayList<>(numBreakpoints);

        for (final Tuple2<Tuple2<String,String>, BreakpointAlignment> assembledBreakpointPair : assembledBreakpointsList) {
            final BreakpointAlignment breakpointAlignment = assembledBreakpointPair._2;
            numAssembledBreakpoints = numAssembledBreakpoints + 1;
            final int assembledBreakpointMapq = Math.min(breakpointAlignment.region1.mapqual, breakpointAlignment.region2.mapqual);
            if (assembledBreakpointMapq == 60) {
                highMqMappings = highMqMappings + 1;
            }
            mqs.add(assembledBreakpointMapq);
            final int assembledBreakpointAlignmentLength =
                    Math.min(breakpointAlignment.region1.referenceInterval.size(),
                            breakpointAlignment.region2.referenceInterval.size()) - breakpointAlignment.region1.overlapOnContig(breakpointAlignment.region2);
            alignLengths.add(assembledBreakpointAlignmentLength);
            maxAlignLength = Math.max(maxAlignLength, assembledBreakpointAlignmentLength);
            breakpointIds.add(assembledBreakpointPair._1._1);
            assembledContigIds.add(breakpointAlignment.contigId);
            insertionMappings.addAll(breakpointAlignment.insertionMappings);
        }

        return createVariant(numAssembledBreakpoints, highMqMappings, mqs, alignLengths, maxAlignLength, breakpointAllele, breakpointIds, assembledContigIds, broadcastReference.getValue(), insertionMappings);
    }

    @VisibleForTesting
    private static VariantContext createVariant(final int numAssembledBreakpoints, final int highMqMappings, final List<Integer> mqs, final List<Integer> alignLengths, final int maxAlignLength,
                                                final BreakpointAllele breakpointAllele, final List<String> breakpointIds, final List<String> assembledContigIds, final ReferenceMultiSource reference,
                                                final List<String> insertionMappings) throws IOException {
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
                .attribute(GATKSVVCFHeaderLines.ASSEMBLY_IDS, breakpointIds.stream().collect(Collectors.joining(",")))
                .attribute(GATKSVVCFHeaderLines.CONTIG_IDS, assembledContigIds.stream().map(s -> s.replace(" ", "_")).collect(Collectors.joining(",")));

        if (breakpointAllele.insertedSequence.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, breakpointAllele.insertedSequence);
        }

        if (insertionMappings.size() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, insertionMappings.stream().collect(Collectors.joining("|")));
        }

        if (breakpointAllele.homology.length() > 0) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY, breakpointAllele.homology);
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.HOMOLOGY_LENGTH, breakpointAllele.homology.length());
        }

        if (breakpointAllele.getInversionType() == INV_5_TO_3) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_5_TO_3, "");
        }

        if (breakpointAllele.getInversionType() == INV_3_TO_5) {
            vcBuilder = vcBuilder.attribute(GATKSVVCFHeaderLines.INV_3_TO_5, "");
        }


        return vcBuilder.make();
    }


    private static String getInversionId(final BreakpointAllele breakpointAllele) {
        final String invType = breakpointAllele.getInversionType().name();
        return invType + "_" + breakpointAllele.leftAlignedLeftBreakpoint.getContig() + "_" + breakpointAllele.leftAlignedLeftBreakpoint.getStart() + "_" + breakpointAllele.leftAlignedRightBreakpoint.getStart();
    }


    static Tuple2<BreakpointAllele, Tuple2<Tuple2<String,String>, BreakpointAlignment>> keyByBreakpointAllele(final Tuple2<Tuple2<String, String>, BreakpointAlignment> breakpointIdAndAssembledBreakpoint) {
        final BreakpointAllele breakpointAllele = breakpointIdAndAssembledBreakpoint._2.getBreakpointAllele();
        return new Tuple2<>(breakpointAllele, new Tuple2<>(breakpointIdAndAssembledBreakpoint._1, breakpointIdAndAssembledBreakpoint._2));
    }

}
