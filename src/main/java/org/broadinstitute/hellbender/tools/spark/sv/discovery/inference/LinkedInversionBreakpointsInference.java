package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.InversionBreakendPreFilter.OverlappingPair;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * Depending on how the two intervals overlap each other (see {@link OverlappingScenario}),
 * the signatures tell us what possible types of events are involved.
 * However, there are possibly some degeneracies in some of the overlap scenarios
 * that needs to be resolved.
 */
public final class LinkedInversionBreakpointsInference implements Serializable {
    private static final long serialVersionUID = 1L;
    private static boolean DEV_MODE = false;

    // see StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.allowedShortFragmentOverhang
    private static final int ALLOWED_SHORT_FRAGMENT_OVERHANG = 10;

    private static final int PADDING_LENGTH = 6 * 151;     // length of short read assumed to be 151, to be used for padding when constructing artificial reference

    private static final int tempRefFastaLineLength = 100;  // line length for writing out temporary fasta files containing artificial reference sequences

    // TODO: 5/16/18 to be moved to GATKSVVCFConstants
    static final String INV_FLANKING_MICRO_DEL_KEY = "FLANKING_MICRO_DEL";
    static final String INV_FLANKING_MICRO_DUP_KEY = "FLANKING_MICRO_DUP";
    static final String INV_SOURCE_BND_EVENT_KEY = "SOURCE_BNDS";
    static final String EVENT_TYPE_KEY = "EVENT_TYPE";

    // TODO: 5/16/18 to be moved to SimpleSVType
    private static final Allele INV_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INV.name()), false);
    private static final Allele DEL_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.DEL.name()), false);
    private static final Allele INS_SYMB_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SimpleSVType.SupportedType.INS.name()), false);

    private static final String DISPERSED_INVERTED_DUPLICATION_UPSTREAM_INSERTION_ORIENTATIONS = "-+";
    private static final String DISPERSED_INVERTED_DUPLICATION_DOWNSTREAM_INSERTION_ORIENTATIONS = "+-";

    private static String makeID(final String typeName, final String chr, final int start, final int stop) {
        return typeName + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + chr + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + start + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + stop;
    }

    public enum OverlappingScenario {
        // what it means for possible variation(s) is described in toString()
        THREE_CONTAINS_FIVE,      // interval from INV33 contains that from INV55
        FIVE_CONTAINS_THREE,      // interval from INV55 contains that from INV33
        THREE_INTERLEAVES_FIVE,   // interval from INV33 upstream of that from INV55, but doesn't contain it
        FIVE_INTERLEAVES_THREE;   // interval from INV55 upstream of that from INV33, but doesn't contain it

        @Override
        public String toString() { // a "primed-annotated" block is inverted
            switch (this) {
                case THREE_CONTAINS_FIVE:
                    return " ABC -> A(B|B')A' ";
                case FIVE_CONTAINS_THREE:
                    return " ABC -> C'(B|B')C ";
                case THREE_INTERLEAVES_FIVE:
                    return " ABC -> AC'B'A'C  ";
                case FIVE_INTERLEAVES_THREE:
                    return " ABC ->    B'     ";
                default: throw new IllegalStateException();
            }
        }

        String getDescription() {
            return toString();
        }
    }

    //==================================================================================================================

    /**
     * Main interface.
     * Handles overlapping pairs of breakends.
     * @see OverlappingScenario
     */
    public static List<VariantContext> makeInterpretation(final List<OverlappingPair> overlappingPairs,
                                                          final String pathToFastqsDir,
                                                          final ReferenceMultiSource reference,
                                                          final Logger toolLogger) {

        // two passes on the pairs, one pass to extract relevant short reads, one pass to make inference
        final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2Reads = extractReads(overlappingPairs, pathToFastqsDir);

        return overlappingPairs.stream()
                .map(LinkedInversionBreakpointsInference::factory)
                .flatMap(overlappingPairHandler ->
                        overlappingPairHandler.toVariants(reference, assemblyID2Reads, pathToFastqsDir).stream())
                .collect(Collectors.toList());
    }

    //==================================================================================================================

    private abstract static class OverlappingPairHandler implements Serializable {
        private static final long serialVersionUID = 1L;

        protected final OverlappingPair overlappingPair;

        OverlappingPairHandler(final OverlappingPair overlappingPair) {
            this.overlappingPair = overlappingPair;
        }

        /**
         * Inference logic to be implemented case by case.
         * @param reference         reference for extracting relevant bases
         * @param assemblyID2Reads  map from assembly ID to short reads used for its assembly process
         * @param pathToFastqsDir   a development artifact, because we want to save alignment results to file and study
         * @return                  interpreted variants.
         */
        abstract List<VariantContext> toVariants(final ReferenceMultiSource reference,
                                                 final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2Reads,
                                                 final String pathToFastqsDir);
    }

    /**
     * Depending on overlapping signatures of input {@code overlappingPair},
     * make appropriate handler using scenarios recognized in {@link OverlappingScenario}.
     */
    private static OverlappingPairHandler factory(final OverlappingPair overlappingPair) {

        final SimpleInterval intervalFive = overlappingPair.fivePrimeBreakEnd.getSpanningInterval();
        final SimpleInterval intervalThree = overlappingPair.threePrimeBreakEnd.getSpanningInterval();

        if (intervalFive.contains(intervalThree)) {
            return new FiveContainsThree(overlappingPair);
        } else if (intervalThree.contains(intervalFive)) {
            return new ThreeContainsFive(overlappingPair);
        } else if (compareIntervals(intervalFive, intervalThree) <= 0) {
            return new FiveInterleavesThree(overlappingPair);
        } else {
            return new ThreeInterleavesFive(overlappingPair);
        }
    }

    //==================================================================================================================

    /**
     * Scenario where re-alignment of short reads back to artificial reference haplotypes
     * is NOT necessary
     * to resolve which variants are here.
     */
    abstract static class NonDegenerateScenario extends OverlappingPairHandler {
        private static final long serialVersionUID = 1L;

        NonDegenerateScenario(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }
    }
    
    private static final class FiveInterleavesThree extends NonDegenerateScenario {
        private static final long serialVersionUID = 1L;
        FiveInterleavesThree(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }

        // TODO: 5/7/18 refactor the following code and learn from the refactoring made in PR #4602
        @Override
        public List<VariantContext> toVariants(final ReferenceMultiSource reference,
                                               final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2Reads,
                                               final String pathToFastqsDir) {

            final List<String> contigNames =
                    Stream.concat(overlappingPair.fivePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES),
                            overlappingPair.threePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES))
                            .sorted().distinct().collect(Collectors.toList());

            final String sourceIDs = overlappingPair.fivePrimeBreakEnd.getID().replaceAll("_1$", "")
                    + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR
                    + overlappingPair.threePrimeBreakEnd.getID().replaceAll("_1$", "");


            final List<VariantContextBuilder> result = new ArrayList<>(3);


            byte[] refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr,
                    overlappingPair.threeIntervalLeftBoundary, overlappingPair.threeIntervalLeftBoundary));
            final SimpleInterval invertedInterval = new SimpleInterval(overlappingPair.chr,
                    overlappingPair.threeIntervalLeftBoundary, overlappingPair.fiveIntervalRightBoundary);
            final VariantContextBuilder inversionBuilder = makeInversion(invertedInterval, Allele.create(refBase, true));

            if ( overlappingPair.threeIntervalLeftBoundary - overlappingPair.fiveIntervalLeftBoundary - 1 > 49 ) {
                refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary, overlappingPair.fiveIntervalLeftBoundary));
                final SimpleInterval deletedInterval = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.fiveIntervalLeftBoundary, overlappingPair.threeIntervalLeftBoundary);
                final VariantContextBuilder leftDel = makeDeletion(deletedInterval, Allele.create(refBase, true));
                result.add(0, leftDel);
            } else if (overlappingPair.fiveIntervalLeftBoundary + 1 != overlappingPair.threeIntervalLeftBoundary) {
                inversionBuilder.attribute(INV_FLANKING_MICRO_DEL_KEY,
                        new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary + 1, overlappingPair.threeIntervalLeftBoundary - 1).toString());
            }

            if ( overlappingPair.threeIntervalRightBoundary - overlappingPair.fiveIntervalRightBoundary - 1 > 49 ) {
                refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary, overlappingPair.fiveIntervalLeftBoundary));
                final SimpleInterval deletedInterval = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.fiveIntervalRightBoundary, overlappingPair.threeIntervalRightBoundary);
                final VariantContextBuilder rightDel = makeDeletion(deletedInterval, Allele.create(refBase, true));
                result.add(rightDel);
            } else if (overlappingPair.fiveIntervalRightBoundary + 1 != overlappingPair.threeIntervalRightBoundary) {
                inversionBuilder.attribute(INV_FLANKING_MICRO_DEL_KEY,
                        new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary + 1, overlappingPair.threeIntervalRightBoundary - 1).toString());
            }

            result.add(inversionBuilder); // insertion order doesn't matter
            return result.stream().map(builder ->
                    builder.attribute(GATKSVVCFConstants.CONTIG_NAMES, contigNames)
                            .attribute(INV_SOURCE_BND_EVENT_KEY, sourceIDs).make())
                    .collect(Collectors.toList());
        }
    }

    private static final class ThreeInterleavesFive extends NonDegenerateScenario {
        private static final long serialVersionUID = 1L;
        ThreeInterleavesFive(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }

        @Override
        public List<VariantContext> toVariants(final ReferenceMultiSource reference,
                                               final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2Reads,
                                               final String pathToFastqsDir) {

            final List<String> contigNames =
                    Stream.concat(overlappingPair.threePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES),
                            overlappingPair.fivePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES))
                            .sorted().distinct().collect(Collectors.toList());

            final String sourceIDs = overlappingPair.threePrimeBreakEnd.getID().replaceAll("_1$", "")
                    + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR
                    + overlappingPair.fivePrimeBreakEnd.getID().replaceAll("_1$", "");

            final List<VariantContextBuilder> result = new ArrayList<>(3);

            byte[] refBase = getRefBase(reference, SVLocalContext.makeOneBpInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary));
            final SimpleInterval invertedInterval = new SimpleInterval(overlappingPair.chr,
                    overlappingPair.fiveIntervalLeftBoundary - 1, overlappingPair.threeIntervalRightBoundary - 1);
            final VariantContextBuilder inversionBuilder = makeInversion(invertedInterval, Allele.create(refBase, true));

            int svLen = overlappingPair.fiveIntervalLeftBoundary - overlappingPair.threeIntervalLeftBoundary + 1;
            SimpleInterval dupRange = new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary, overlappingPair.fiveIntervalLeftBoundary);
            if ( svLen > 49 ) {

                refBase = getRefBase(reference, SVLocalContext.makeOneBpInterval(overlappingPair.chr, overlappingPair.threeIntervalRightBoundary - 1));
                final VariantContextBuilder leftDup = makeDispersedDuplication(overlappingPair.chr,
                        overlappingPair.threeIntervalRightBoundary - 1,
                        dupRange,DISPERSED_INVERTED_DUPLICATION_DOWNSTREAM_INSERTION_ORIENTATIONS,
                        Allele.create(refBase, true));
                result.add(leftDup);
            } else if (overlappingPair.threeIntervalLeftBoundary + 1 != overlappingPair.fiveIntervalLeftBoundary) {
                inversionBuilder.attribute(INV_FLANKING_MICRO_DUP_KEY, dupRange.toString());
            }

            svLen = overlappingPair.fiveIntervalRightBoundary - overlappingPair.threeIntervalRightBoundary + 1;
            dupRange = new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalRightBoundary, overlappingPair.fiveIntervalRightBoundary);
            if ( svLen > 49 ) {
                refBase = getRefBase(reference, SVLocalContext.makeOneBpInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary));
                final VariantContextBuilder rightDup = makeDispersedDuplication(overlappingPair.chr,
                        overlappingPair.fiveIntervalLeftBoundary,
                        dupRange, DISPERSED_INVERTED_DUPLICATION_UPSTREAM_INSERTION_ORIENTATIONS,
                        Allele.create(refBase, true));
                result.add(rightDup);
            } else if (overlappingPair.threeIntervalRightBoundary + 1 != overlappingPair.fiveIntervalRightBoundary) {
                inversionBuilder.attribute(INV_FLANKING_MICRO_DUP_KEY, dupRange.toString());
            }

            result.add(inversionBuilder); // insertion order doesn't matter
            return result.stream().map(builder ->
                    builder.attribute(GATKSVVCFConstants.CONTIG_NAMES, contigNames)
                            .attribute(INV_SOURCE_BND_EVENT_KEY, sourceIDs).make())
                    .collect(Collectors.toList());
        }
    }

    //==================================================================================================================

    /**
     * Scenario where re-alignment of short reads back to artificial reference haplotypes
     * is necessary
     * to resolve which variants are here.
     */
    abstract static class DegenerateScenario extends OverlappingPairHandler {
        private static final long serialVersionUID = 1L;
        DegenerateScenario(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }

        @Override
        public final List<VariantContext> toVariants(final ReferenceMultiSource reference,
                                                     final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2ShortReads,
                                                     final String pathToFastqsDir) {

            final Set<String> supportingAssemblyIDs = new HashSet<>( overlappingPair.fivePrimeBreakEnd.getSupportingAssemblyIDs() );
            supportingAssemblyIDs.addAll(overlappingPair.threePrimeBreakEnd.getSupportingAssemblyIDs());

            final Tuple2<AlignmentResult.ArtificialReferenceHaplotypes, AlignmentResult.ArtificialReferenceHaplotypes>
                    artificialReferenceGroups = makeTwoArtificialReferenceGroups(reference);

            final List<AlignmentResult> alignmentsForGroupOne =
                    alignShortReadsToArtificialReference(artificialReferenceGroups._1, assemblyID2ShortReads, supportingAssemblyIDs,
                            DEV_MODE ? makeReferenceFileName(pathToFastqsDir, supportingAssemblyIDs, true) : null);
            final List<AlignmentResult> alignmentsForGroupTwo =
                    alignShortReadsToArtificialReference(artificialReferenceGroups._2, assemblyID2ShortReads, supportingAssemblyIDs,
                            DEV_MODE ? makeReferenceFileName(pathToFastqsDir, supportingAssemblyIDs, false) : null);

            if (DEV_MODE) {
                writeAlignments(pathToFastqsDir, alignmentsForGroupOne, true);
                writeAlignments(pathToFastqsDir, alignmentsForGroupTwo, false);
            }
            return breakDegeneracyAndInterpret(alignmentsForGroupOne, alignmentsForGroupTwo, reference);
        }

        private String makeReferenceFileName(final String workingDir,
                                             final Set<String> supportingAssemblyIDs,
                                             final boolean simpler) {

            String postfix = String.join("_", supportingAssemblyIDs);

            return workingDir + (workingDir.endsWith("/") ? "" : "/") + postfix + (simpler ? "_0" : "_1") + ".fasta";
        }

        protected final List<AlignmentResult> alignShortReadsToArtificialReference(final AlignmentResult.ArtificialReferenceHaplotypes artificialReferenceHaplotypes,
                                                                                   final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2ShortReads,
                                                                                   final Set<String> supportingAssemblyIDs,
                                                                                   final String pathToSaveArtificialReferenceFasta) {

            try {
                final List<AlignmentResult> alignmentResults = new ArrayList<>();

                try (final ReferenceSequencesAligner aligner =
                             new ReferenceSequencesAligner(
                                     artificialReferenceHaplotypes.getRefNames(),
                                     artificialReferenceHaplotypes.getDescribedRefContigs(),
                                     true, true, tempRefFastaLineLength,
                                     pathToSaveArtificialReferenceFasta) ) {

                    final SAMFileHeader samFileHeader = new SAMFileHeader(aligner.getDict());

                    final Map<String, List<List<SAMRecord>>> assemblyID2Sam = new HashMap<>(supportingAssemblyIDs.size());
                    for (final String assemblyID : supportingAssemblyIDs) {
                        final List<SVFastqUtils.FastqRead> fastqReads = assemblyID2ShortReads.get(assemblyID);
                        final List<List<SAMRecord>> samRecords = aligner.align(fastqReads, samFileHeader, true);
                        assemblyID2Sam.put(assemblyID, samRecords);
                    }

                    alignmentResults.add( new AlignmentResult(artificialReferenceHaplotypes, samFileHeader, assemblyID2Sam) );
                }
                return alignmentResults;
            }
            catch (final IOException ioex) {
                throw new GATKException("fdsfs", ioex);
            }
        }
        
        abstract Tuple2<AlignmentResult.ArtificialReferenceHaplotypes, AlignmentResult.ArtificialReferenceHaplotypes>
        makeTwoArtificialReferenceGroups(final ReferenceMultiSource reference);

        abstract List<VariantContext> breakDegeneracyAndInterpret(final List<AlignmentResult> alignmentsForGroupOne,
                                                                  final List<AlignmentResult> alignmentsForGroupTwo,
                                                                  final ReferenceMultiSource reference);

        private static void writeAlignments(final String workingDir, final List<AlignmentResult> alignmentResults,
                                            final boolean simpler) {

            for (final AlignmentResult alignmentResult : alignmentResults) {

                final SAMFileHeader samFileHeader = alignmentResult.samFileHeader;
                SAMFileHeader clone = samFileHeader.clone();
                clone.setSortOrder(SAMFileHeader.SortOrder.coordinate);
                final SAMRecordComparator localComparator = new SAMRecordCoordinateComparator();

                for (final Map.Entry<String, List<List<SAMRecord>>> pair : alignmentResult.assemblyID2Sam.entrySet()) {
                    String key = pair.getKey() + (simpler ? "_0" : "_1");
                    final String bamOut = workingDir.endsWith("/") ? workingDir.concat(key).concat(".bam")
                                                                   : workingDir.concat("/").concat(key).concat(".bam");

                    final List<SAMRecord> samRecords = pair.getValue().stream().flatMap(List::stream).collect(Collectors.toList());
                    samRecords.sort(localComparator);
                    SVFileUtils.writeSAMFile(bamOut, samRecords.iterator(), clone, true);
                }
            }
        }

        private enum PairOrientation {
            LL, RR, RL,
            LR
        }
        static final PairOrientation getPairOrientation(final SAMRecord samRecord) {
            if (samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
                return samRecord.getReadNegativeStrandFlag() ? PairOrientation.LL : PairOrientation.RR;
            else if (samRecord.getReadNegativeStrandFlag() ?
                    samRecord.getStart() + ALLOWED_SHORT_FRAGMENT_OVERHANG < samRecord.getMateAlignmentStart() :
                    samRecord.getStart() - ALLOWED_SHORT_FRAGMENT_OVERHANG > samRecord.getMateAlignmentStart()) {
                return PairOrientation.RL;
            } else {
                return PairOrientation.LR;
            }
        }

        final Long[] getLeftyAndRightPairCounts(final List<AlignmentResult> alignmentsForGroupOne,
                                                final List<AlignmentResult> alignmentsForGroupTwo) {
            final EnumMap<PairOrientation, Long> simple = alignmentsForGroupOne.stream()
                    .flatMap(alignmentResult -> alignmentResult.assemblyID2Sam.values().stream().flatMap(List::stream).flatMap(List::stream))
                    .filter(samRecord -> !(samRecord.getReadUnmappedFlag() || samRecord.isSecondaryOrSupplementary() || samRecord.getMateUnmappedFlag()))
                    .collect(Collectors.groupingBy(DegenerateScenario::getPairOrientation,
                            () -> new EnumMap<>(PairOrientation.class),
                            Collectors.counting()));

            final EnumMap<PairOrientation, Long> complex = alignmentsForGroupTwo.stream()
                    .flatMap(alignmentResult -> alignmentResult.assemblyID2Sam.values().stream().flatMap(List::stream).flatMap(List::stream))
                    .filter(samRecord -> !(samRecord.getReadUnmappedFlag() || samRecord.isSecondaryOrSupplementary() || samRecord.getMateUnmappedFlag()))
                    .collect(Collectors.groupingBy(DegenerateScenario::getPairOrientation,
                            () -> new EnumMap<>(PairOrientation.class),
                            Collectors.counting()));

            final String assemblyIDs = String.join(",",
                    alignmentsForGroupOne.stream().flatMap(ar -> ar.getAssemblyIDs().stream()).collect(Collectors.toSet()));

            final Long simpleCaseLefties = simple.get(PairOrientation.LL) == null ? 0L : simple.get(PairOrientation.LL);
            final Long simpleCaseRighties = simple.get(PairOrientation.RR) == null ? 0L : simple.get(PairOrientation.RR);
            final Long complexCaseLefties = complex.get(PairOrientation.LL) == null ? 0L : complex.get(PairOrientation.LL);
            final Long complexCaseRighties = complex.get(PairOrientation.RR) == null ? 0L : complex.get(PairOrientation.RR);
            if (DEV_MODE) {
                String line = "BED" + "\t" + overlappingPair.chr + "\t" +
                        Math.min(overlappingPair.fiveIntervalLeftBoundary, overlappingPair.threeIntervalLeftBoundary) + "\t" +
                        Math.max(overlappingPair.fiveIntervalRightBoundary, overlappingPair.threeIntervalRightBoundary) + "\t" +
                        simpleCaseLefties + "," + simpleCaseRighties + ";" + complexCaseLefties + "," + complexCaseRighties + ";" +
                        assemblyIDs;
                System.out.println(line);
            }
            return new Long[]{simpleCaseLefties, simpleCaseRighties, complexCaseLefties, complexCaseRighties};
        }

        // TODO: 5/18/18 couldn't be more naive, but is working
        boolean shouldTriggerInversionCall(final List<AlignmentResult> alignmentsForGroupOne,
                                           final List<AlignmentResult> alignmentsForGroupTwo) {
            final Long[] leftyAndRightPairCounts = getLeftyAndRightPairCounts(alignmentsForGroupOne, alignmentsForGroupTwo);
            return leftyAndRightPairCounts[0] + leftyAndRightPairCounts[1] >= 30;
        }
    }

    /**
     * Holding results of aligning short reads to artificial references.
     */
    private static final class AlignmentResult {
        final ArtificialReferenceHaplotypes artificialReferenceHaplotypes;
        final SAMFileHeader samFileHeader;
        final Map<String, List<List<SAMRecord>>> assemblyID2Sam;

        AlignmentResult(final ArtificialReferenceHaplotypes artificialReferenceHaplotypes,
                        final SAMFileHeader samFileHeader,
                        final Map<String, List<List<SAMRecord>>> alignments) {
            this.artificialReferenceHaplotypes = artificialReferenceHaplotypes;
            this.samFileHeader = samFileHeader;
            this.assemblyID2Sam = alignments;
        }

        private Set<String> getAssemblyIDs() {
            return assemblyID2Sam.keySet();
        }

        private static class ArtificialReferenceHaplotypes {

            private final List<ReferenceSequencesAligner.DescribedRefContig> describedRefContigs;

            private ArtificialReferenceHaplotypes(final List<ReferenceSequencesAligner.DescribedRefContig> describedRefContigs) {
                this.describedRefContigs = describedRefContigs;
            }

            List<String> getRefNames() {
                return describedRefContigs.stream().map(ReferenceSequencesAligner.DescribedRefContig::getName).collect(Collectors.toList());
            }

            List<ReferenceSequencesAligner.DescribedRefContig> getDescribedRefContigs() {
                return describedRefContigs;
            }
        }
    }

    /**
     * As described in {@link OverlappingScenario#toString()},
     * there's a degeneracy in terms of whether the B block is inverted.
     * Since C'BC and C'B'C are reverse complement of each other,
     * we need to check alignment signatures from short reads
     * that link this middle group of to left/right flanking regions
     * to break this degeneracy.
     */
    private static final class ThreeContainsFive extends DegenerateScenario {
        private static final long serialVersionUID = 1L;
        ThreeContainsFive(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }

        /**
         * Based on the pair interval overlap scenario,
         * artificial reference haplotypes are made for
         * breaking degeneracy.
         * (again, see {@link OverlappingScenario#toString()})
         */
        @Override
        Tuple2<AlignmentResult.ArtificialReferenceHaplotypes, AlignmentResult.ArtificialReferenceHaplotypes>
        makeTwoArtificialReferenceGroups(final ReferenceMultiSource reference) {

            try {
                final ReferenceBases leftFlanking = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary - PADDING_LENGTH,
                                                                                    overlappingPair.threeIntervalLeftBoundary));
                final ReferenceBases rightFlankingRegion = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalRightBoundary,
                                                                         overlappingPair.threeIntervalRightBoundary + PADDING_LENGTH));
                final ReferenceBases blockB = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary, overlappingPair.fiveIntervalRightBoundary));
                final byte[] invertedBlockB = Arrays.copyOf(blockB.getBases(), blockB.getBases().length);
                SequenceUtil.reverseComplement(invertedBlockB);

//                final boolean blockCDeleted = overlappingPair.threeIntervalRightBoundary - overlappingPair.fiveIntervalRightBoundary > 1;
                // if blockCDeleted is true, then then the deleted block C would be calculated as follows:
//                ReferenceBases blockC = reference
//                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary + 1, overlappingPair.threeIntervalRightBoundary - 1));

                try (final ByteArrayOutputStream outputStream = new ByteArrayOutputStream()) {
                    // guard against edge case
                    final boolean blockAInverseDisperseDuplicated = (overlappingPair.fiveIntervalLeftBoundary - overlappingPair.threeIntervalLeftBoundary > 1);

                    outputStream.write(leftFlanking.getBases());

                    final ReferenceBases blockA;
                    final byte[] invertedBlockA;
                    String description;
                    String name;
                    if (blockAInverseDisperseDuplicated) {
                        blockA = reference
                                .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary - 1, overlappingPair.fiveIntervalLeftBoundary - 1));
                        invertedBlockA = Arrays.copyOf(blockA.getBases(), blockA.getBases().length);
                        SequenceUtil.reverseComplement(invertedBlockA);

                        outputStream.write(blockA.getBases());
                        outputStream.write(blockB.getBases());
                        outputStream.write(invertedBlockA);

                        description = leftFlanking.getInterval() + "+" + blockA.getInterval() + "+" +
                                blockB.getInterval() + "+(" + blockA.getInterval() + ")'+" + rightFlankingRegion.getInterval();
                        name = "LABaR";
                    } else {
                        blockA = null;
                        invertedBlockA = null;
                        outputStream.write(blockB.getBases());

                        description = leftFlanking.getInterval() + "+" +
                                blockB.getInterval() + "+" + rightFlankingRegion.getInterval();
                        name = "LBR";
                    }

                    outputStream.write(rightFlankingRegion.getBases());
                    final byte[] sequencesForFirstGroup = outputStream.toByteArray();
                    outputStream.flush();

                    final ReferenceSequencesAligner.DescribedRefContig firstGroup =
                            new ReferenceSequencesAligner.DescribedRefContig(
                                    name, description, sequencesForFirstGroup
                            );


                    outputStream.reset();
                    outputStream.write(leftFlanking.getBases());

                    if (blockAInverseDisperseDuplicated) {
                        outputStream.write(blockA.getBases());
                        outputStream.write(invertedBlockB);
                        outputStream.write(invertedBlockA);
                        description = leftFlanking.getInterval() + "+" + blockA.getInterval() + "+(" +
                                blockB.getInterval() + ")'+(" + blockA.getInterval() + ")'+" + rightFlankingRegion.getInterval();
                        name = "LAbaR";
                    } else {
                        outputStream.write(invertedBlockB);
                        description = leftFlanking.getInterval() + "+(" +
                                blockB.getInterval() + ")'+" + rightFlankingRegion.getInterval();
                        name = "LbR";
                    }
                    outputStream.write(rightFlankingRegion.getBases());
                    final byte[] sequencesForSecondGroup = outputStream.toByteArray();
                    outputStream.flush();

                    final ReferenceSequencesAligner.DescribedRefContig secondGroup = new ReferenceSequencesAligner.DescribedRefContig(
                            name, description, sequencesForSecondGroup
                    );

                    return new Tuple2<>(new AlignmentResult.ArtificialReferenceHaplotypes(Collections.singletonList(firstGroup)),
                                        new AlignmentResult.ArtificialReferenceHaplotypes(Collections.singletonList(secondGroup)));
                }
            } catch (final IOException ioex) {
                throw new GATKException("Can not get reference bases", ioex);
            }
        }

        @Override
        List<VariantContext> breakDegeneracyAndInterpret(final List<AlignmentResult> alignmentsForGroupOne,
                                                         final List<AlignmentResult> alignmentsForGroupTwo,
                                                         final ReferenceMultiSource reference) {
            final List<String> contigNames =
                    Stream.concat(overlappingPair.fivePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES),
                            overlappingPair.threePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES))
                            .sorted().distinct().collect(Collectors.toList());

            final String sourceIDs = overlappingPair.fivePrimeBreakEnd.getID().replaceAll("_1$", "")
                    + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR
                    + overlappingPair.threePrimeBreakEnd.getID().replaceAll("_1$", "");

            final List<VariantContextBuilder> result = new ArrayList<>(3);
            if (overlappingPair.threeIntervalRightBoundary - overlappingPair.fiveIntervalRightBoundary - 1 > 49) {
                final SimpleInterval delRange = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.fiveIntervalRightBoundary, overlappingPair.threeIntervalRightBoundary - 1);
                byte[] refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary,
                        overlappingPair.fiveIntervalRightBoundary));
                final VariantContextBuilder builder = makeDeletion(delRange, Allele.create(refBase, true));
                result.add(builder);
            }
            if (overlappingPair.fiveIntervalLeftBoundary - overlappingPair.threeIntervalLeftBoundary + 1 > 49) {
                final SimpleInterval dupRange = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.threeIntervalLeftBoundary, overlappingPair.fiveIntervalLeftBoundary);
                byte[] refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary,
                        overlappingPair.fiveIntervalRightBoundary));
                final VariantContextBuilder builder = makeDispersedDuplication(
                        overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary,
                        dupRange, DISPERSED_INVERTED_DUPLICATION_DOWNSTREAM_INSERTION_ORIENTATIONS,
                        Allele.create(refBase, true));
                result.add(builder);
            }

            if (shouldTriggerInversionCall(alignmentsForGroupOne, alignmentsForGroupTwo)) {
                final SimpleInterval invRange = new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary + 1,
                        overlappingPair.fiveIntervalRightBoundary);
                final VariantContextBuilder builder = makeInversion(invRange, overlappingPair.fivePrimeBreakEnd.getReference());
                result.add(builder);
            }

            return result.stream()
                    .map(builder ->
                            builder.attribute(GATKSVVCFConstants.CONTIG_NAMES, contigNames)
                                    .attribute(INV_SOURCE_BND_EVENT_KEY, sourceIDs).make())
                    .collect(Collectors.toList());
        }
    }

    /**
     * As described in {@link OverlappingScenario#toString()},
     * there's a degeneracy in terms of whether the B block is inverted.
     * Since ABA' and AB'A' are reverse complement of each other,
     * we need to check alignment signatures from short reads
     * that link this middle group of to left/right flanking regions
     * to break this degeneracy.
     */
    private static final class FiveContainsThree extends DegenerateScenario {
        private static final long serialVersionUID = 1L;
        FiveContainsThree(final OverlappingPair overlappingPair) {
            super(overlappingPair);
        }

        /**
         * Based on the pair interval overlap scenario,
         * artificial reference haplotypes are made for
         * breaking degeneracy.
         * (again, see {@link OverlappingScenario#toString()})
         */
        @Override
        Tuple2<AlignmentResult.ArtificialReferenceHaplotypes, AlignmentResult.ArtificialReferenceHaplotypes>
        makeTwoArtificialReferenceGroups(final ReferenceMultiSource reference) {

            try {
                final ReferenceBases leftFlanking = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary - PADDING_LENGTH,
                                overlappingPair.fiveIntervalLeftBoundary));
                final ReferenceBases rightFlankingRegion = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalRightBoundary,
                                overlappingPair.fiveIntervalRightBoundary + PADDING_LENGTH));

//                final boolean blockADeleted = overlappingPair.threeIntervalLeftBoundary - overlappingPair.fiveIntervalLeftBoundary > 1;
                // if blockADeleted is true, then then the deleted block A would be calculated as follows:
//                ReferenceBases blockA = reference
//                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary + 1, overlappingPair.threeIntervalLeftBoundary - 1));

                final ReferenceBases blockB = reference
                        .getReferenceBases(new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary, overlappingPair.threeIntervalRightBoundary));
                final byte[] invertedBlockB = Arrays.copyOf(blockB.getBases(), blockB.getBases().length);
                SequenceUtil.reverseComplement(invertedBlockB);

                try (final ByteArrayOutputStream outputStream = new ByteArrayOutputStream() ){
                    // guard against edge case
                    final boolean blockCInverseDisperseDuplicated = (overlappingPair.fiveIntervalRightBoundary - overlappingPair.threeIntervalRightBoundary > 1);

                    outputStream.write(leftFlanking.getBases());

                    final ReferenceBases blockC;
                    final byte[] invertedBlockC;
                    String name;
                    String description;
                    if (blockCInverseDisperseDuplicated) {
                        blockC = reference.getReferenceBases(new SimpleInterval(overlappingPair.chr,
                                overlappingPair.threeIntervalRightBoundary + 1,
                                overlappingPair.fiveIntervalRightBoundary - 1));
                        invertedBlockC = Arrays.copyOf(blockC.getBases(), blockC.getBases().length);
                        SequenceUtil.reverseComplement(invertedBlockC);
                        outputStream.write(invertedBlockC);
                        outputStream.write(blockB.getBases());
                        outputStream.write(blockC.getBases());

                        name = "LcBCR";
                        description = leftFlanking.getInterval() + "+(" + blockC.getInterval() + ")'+" +
                                blockB.getInterval() + "+" + blockC.getInterval()+ "+" +rightFlankingRegion.getInterval();
                    } else {
                        blockC = null;
                        invertedBlockC = null;
                        outputStream.write(blockB.getBases());

                        name = "LBR";
                        description = leftFlanking.getInterval() + "+" +
                                blockB.getInterval() + "+" +rightFlankingRegion.getInterval();
                    }

                    outputStream.write(rightFlankingRegion.getBases());
                    final byte[] sequencesForFirstGroup = outputStream.toByteArray();
                    outputStream.flush();

                    final ReferenceSequencesAligner.DescribedRefContig firstGroup =
                            new ReferenceSequencesAligner.DescribedRefContig(
                                    name, description, sequencesForFirstGroup
                            );

                    outputStream.reset();
                    outputStream.write(leftFlanking.getBases());

                    if (blockCInverseDisperseDuplicated) {
                        outputStream.write(invertedBlockC);
                        outputStream.write(invertedBlockB);
                        outputStream.write(blockC.getBases());

                        name = "LcbCR";
                        description = leftFlanking.getInterval() + "+(" + blockC.getInterval() + ")'+(" +
                                blockB.getInterval() + ")'+" + blockC.getInterval()+ "+" +rightFlankingRegion.getInterval();
                    } else {
                        outputStream.write(invertedBlockB);

                        name = "LbR";
                        description = leftFlanking.getInterval() + "+(" +
                                blockB.getInterval() + ")'+" +rightFlankingRegion.getInterval();
                    }

                    outputStream.write(rightFlankingRegion.getBases());
                    final byte[] sequencesForSecondGroup = outputStream.toByteArray();
                    outputStream.flush();

                    final ReferenceSequencesAligner.DescribedRefContig secondGroup = new ReferenceSequencesAligner.DescribedRefContig(
                            name, description, sequencesForSecondGroup
                    );

                    return new Tuple2<>(new AlignmentResult.ArtificialReferenceHaplotypes(Collections.singletonList(firstGroup)),
                                        new AlignmentResult.ArtificialReferenceHaplotypes(Collections.singletonList(secondGroup)));
                }
            } catch (final IOException ioex) {
                throw new GATKException("Can not get reference bases", ioex);
            }
        }

        @Override
        List<VariantContext> breakDegeneracyAndInterpret(final List<AlignmentResult> alignmentsForGroupOne,
                                                         final List<AlignmentResult> alignmentsForGroupTwo,
                                                         final ReferenceMultiSource reference) {

            final List<String> contigNames =
                    Stream.concat(overlappingPair.fivePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES),
                            overlappingPair.threePrimeBreakEnd.makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES))
                            .sorted().distinct().collect(Collectors.toList());

            final String sourceIDs = overlappingPair.fivePrimeBreakEnd.getID().replaceAll("_1$", "")
                    + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR
                    + overlappingPair.threePrimeBreakEnd.getID().replaceAll("_1$", "");

            final List<VariantContextBuilder> result = new ArrayList<>(3);
            if (overlappingPair.threeIntervalLeftBoundary - overlappingPair.fiveIntervalLeftBoundary - 1 > 49) {
                final SimpleInterval delRange = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.fiveIntervalLeftBoundary, overlappingPair.threeIntervalLeftBoundary - 1);
                final VariantContextBuilder builder = makeDeletion(delRange, overlappingPair.fivePrimeBreakEnd.getReference());
                result.add(builder);
            }
            if (overlappingPair.fiveIntervalRightBoundary - overlappingPair.threeIntervalRightBoundary + 1 > 49) {
                final SimpleInterval dupRange = new SimpleInterval(overlappingPair.chr,
                        overlappingPair.threeIntervalRightBoundary, overlappingPair.fiveIntervalRightBoundary);
                final VariantContextBuilder builder = makeDispersedDuplication(
                        overlappingPair.chr, overlappingPair.fiveIntervalLeftBoundary,
                        dupRange, DISPERSED_INVERTED_DUPLICATION_UPSTREAM_INSERTION_ORIENTATIONS,
                        overlappingPair.fivePrimeBreakEnd.getReference());
                result.add(builder);
            }

            if (shouldTriggerInversionCall(alignmentsForGroupOne, alignmentsForGroupTwo)) {
                final SimpleInterval invRange = new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary, overlappingPair.threeIntervalRightBoundary);
                byte[] refBase = getRefBase(reference, new SimpleInterval(overlappingPair.chr, overlappingPair.threeIntervalLeftBoundary - 1,
                        overlappingPair.threeIntervalLeftBoundary - 1));
                final VariantContextBuilder builder = makeInversion(invRange, Allele.create(refBase, true));
                result.add(builder);
            }

            return result.stream()
                    .map(builder ->
                            builder.attribute(GATKSVVCFConstants.CONTIG_NAMES, contigNames)
                                    .attribute(INV_SOURCE_BND_EVENT_KEY, sourceIDs).make())
                    .collect(Collectors.toList());
        }

    }

    //==================================================================================================================

    private static int compareIntervals(final SimpleInterval first, final SimpleInterval second) {
        // compare start position
        int result = Integer.compare(first.getStart(), second.getStart());
        if (result == 0) {
            // compare end position
            result = Integer.compare(first.getEnd(), second.getEnd());
        }
        return result;
    }

    /**
     * Get short reads used for local assembly that triggered the BND calls in {@code overlappingPairs}.
     */
    private static Map<String, List<SVFastqUtils.FastqRead>> extractReads(final List<OverlappingPair> overlappingPairs,
                                                                          final String pathToFastqsDir) {
        final Set<String> interestingAssemblyIDs = new HashSet<>();
        overlappingPairs.forEach(pair -> {
            interestingAssemblyIDs.addAll( pair.fivePrimeBreakEnd.getSupportingAssemblyIDs() );
            interestingAssemblyIDs.addAll( pair.threePrimeBreakEnd.getSupportingAssemblyIDs() );
        });

        final Map<String, List<SVFastqUtils.FastqRead>> assemblyID2ShortReads = new HashMap<>();
        final String dir = pathToFastqsDir.endsWith("/") ? pathToFastqsDir : pathToFastqsDir + "/";
        try (Stream<Path> paths = Files.walk(IOUtils.getPath(dir))) {
            paths.filter(Files::isRegularFile)
                    .forEach(path ->  {
                        final String fileName = path.toAbsolutePath().toString();
                        int idx = fileName.indexOf(".fastq");
                        if (idx == -1) return; // no a fastq
                        String assemblyId = fileName.substring(idx - 9, idx);// -9 to count chars asm[0-9]{6,6}
                        if (interestingAssemblyIDs.contains(assemblyId)) {
                            String filePath = String.format("%s%s.fastq", dir, assemblyId);
                            List<SVFastqUtils.FastqRead> fastqReads = SVFastqUtils.readFastqFile(filePath);
                            assemblyID2ShortReads.put(assemblyId, fastqReads);
                        }
                    });
        } catch (final IOException ioex) {
            throw new GATKException("Failed to traverse provided directory supposedly containing FASTQ files: " +
                    pathToFastqsDir, ioex);
        }
        return assemblyID2ShortReads;
    }

    private static byte[] getRefBase(final ReferenceMultiSource reference, final SimpleInterval oneBpPos) {
        try {
            return reference.getReferenceBases(oneBpPos).getBases();
        } catch (final IOException ioex) {
            throw new GATKException("Could not read reference for extracting bases", ioex);
        }
    }

    /**
     * Note that {@code delRange} is expected to be pre-process to VCF spec compatible,
     * e.g. if chr1:101-200 is deleted, then {@code delRange} should be chr1:100-200
     */
    @VisibleForTesting
    static VariantContextBuilder makeDeletion(final SimpleInterval delRange, final Allele refAllele) {

        return new VariantContextBuilder()
                .chr(delRange.getContig()).start(delRange.getStart()).stop(delRange.getEnd())
                .alleles(Arrays.asList(refAllele, DEL_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.DEL.name(), delRange.getContig(), delRange.getStart(), delRange.getEnd()))
                .attribute(VCFConstants.END_KEY, delRange.getEnd())
                .attribute(SVLEN, - delRange.size() + 1)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.name());
    }

    @VisibleForTesting
    static VariantContextBuilder makeDispersedDuplication(final String insertionChr, final int insertionPos,
                                                          final SimpleInterval duplicatedRegion, final String duplicationOrientations,
                                                          final Allele refAllele) {
        return new VariantContextBuilder().chr(insertionChr).start(insertionPos).stop(insertionPos)
                .alleles(Arrays.asList(refAllele, INS_SYMB_ALLELE))
                .id("INS-DUPLICATION-DISPERSED" + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + insertionChr + INTERVAL_VARIANT_ID_FIELD_SEPARATOR + insertionPos)
                .attribute(VCFConstants.END_KEY, insertionPos)
                .attribute(SVLEN, duplicatedRegion.size())
                .attribute(SVTYPE, SimpleSVType.SupportedType.INS.name())
                .attribute(EVENT_TYPE_KEY, "DUP:DISPERSED")
                .attribute(DUP_REPEAT_UNIT_REF_SPAN, duplicatedRegion)
                .attribute(DUP_ORIENTATIONS, duplicationOrientations);
    }

    @VisibleForTesting
    static VariantContextBuilder makeInversion(final SimpleInterval invertedRegion, final Allele refAllele) {
        return new VariantContextBuilder()
                .chr(invertedRegion.getContig()).start(invertedRegion.getStart() - 1).stop(invertedRegion.getEnd())     // TODO: 5/2/18 VCF spec doesn't requst left shift by 1 for inversion POS
                .alleles(Arrays.asList(refAllele, INV_SYMB_ALLELE))
                .id(makeID(SimpleSVType.SupportedType.INV.name(), invertedRegion.getContig(), invertedRegion.getStart() - 1, invertedRegion.getEnd()))
                .attribute(VCFConstants.END_KEY, invertedRegion.getEnd())
                .attribute(SVLEN, 0)                                                                 // TODO: 5/2/18 this is following VCF spec,
                .attribute(SVTYPE, SimpleSVType.SupportedType.INV.name());
    }

}
