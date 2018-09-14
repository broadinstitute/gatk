package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputMetaData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;
import scala.Tuple3;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * For extracting simple variants from input GATK-SV complex variants.
 *
 * Some explanation on several concepts:
 *
 * <p>
 *     Anchor ref base:
 *     anchor base is defined per-VCF spec (see 1.4.1#REF version 4.2), that is, for DEL and INS variants
 *     the reference base at the position pointed to by POS, basically:
 *     for DEL, the reference bases immediately following POS are deleted (up to and including the END base),
 *     for INS, the sequence annotated in INSSEQ are inserted immediately after POS.
 * </p>
 *
 * <p>
 *     "Fat" insertion:
 *     they exist because sometimes we have micro deletions surrounding the insertion breakpoint,
 *     so here the strategy is to report them as "fat", i.e. the anchor base and deleted bases are reported in REF;
 *     they are fat in the sense that compared to simple insertions where a single anchor ref base is necessary
 * </p>
 *
 * <p>
 *     It is also assumed that POS and END of the input complex {@link VariantContext} are the boundaries
 *     of the bases where REF and ALT allele share similarity, in other words,
 *     immediately after POS and before END is where the REF and ALT allele differ, and the two path merges at POS/END.
 * </p>
 */
public abstract class SegmentedCpxVariantSimpleVariantExtractor implements Serializable {
    private static final long serialVersionUID = 1L;

    private static int EVENT_SIZE_THRESHOLD = STRUCTURAL_VARIANT_SIZE_LOWER_BOUND - 1;
    private static final String CPX_DERIVED_POSTFIX_STRING = "CPX_DERIVED";
    private static String makeID(final String typeName, final String chr, final int start, final int stop) {
        return typeName + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + chr + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + start + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + stop + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                CPX_DERIVED_POSTFIX_STRING;
    }

    public static final class ExtractedSimpleVariants {
        private final List<VariantContext> reInterpretZeroOrOneSegmentCalls;
        private final List<VariantContext> reInterpretMultiSegmentsCalls;

        private ExtractedSimpleVariants(final List<VariantContext> reInterpretZeroOrOneSegmentCalls,
                                        final List<VariantContext> reInterpretMultiSegmentsCalls) {
            this.reInterpretZeroOrOneSegmentCalls = reInterpretZeroOrOneSegmentCalls;
            this.reInterpretMultiSegmentsCalls = reInterpretMultiSegmentsCalls;
        }

        public List<VariantContext> getReInterpretZeroOrOneSegmentCalls() {
            return reInterpretZeroOrOneSegmentCalls;
        }

        public List<VariantContext> getReInterpretMultiSegmentsCalls() {
            return reInterpretMultiSegmentsCalls;
        }

        public List<VariantContext> getMergedReinterpretedCalls() {
            final ArrayList<VariantContext> merged = new ArrayList<>(reInterpretZeroOrOneSegmentCalls);
            merged.addAll(reInterpretMultiSegmentsCalls);
            return merged;
        }
    }

    // main interface to user code
    public static ExtractedSimpleVariants extract(final JavaRDD<VariantContext> complexVariants,
                                                  final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                                  final JavaRDD<GATKRead> assemblyRawAlignments) {

        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceBroadcast();

        // still does an in-efficient 2-pass on the input RDD: 1 pass for zero- and one-segment calls, the other for multi-segment calls
        // that was due to restriction from how multi-segment calls are to be re-interpreted
        final ZeroAndOneSegmentCpxVariantExtractor zeroAndOneSegmentCpxVariantExtractor = new ZeroAndOneSegmentCpxVariantExtractor();
        final JavaRDD<VariantContext> zeroOrOneSegmentComplexVariants = complexVariants
                .filter(vc -> SVUtils.getAttributeAsStringList(vc, CPX_SV_REF_SEGMENTS).size() < 2)
                .cache();
        final List<VariantContext> reInterpretedZeroAndOneSegmentCalls =
                zeroOrOneSegmentComplexVariants
                        .flatMap(vc -> zeroAndOneSegmentCpxVariantExtractor.extract(vc, referenceBroadcast.getValue()).iterator())
                        .collect();
        zeroOrOneSegmentComplexVariants.unpersist(false);

        final JavaRDD<VariantContext> multiSegmentCalls =
                complexVariants.filter(vc -> SVUtils.getAttributeAsStringList(vc, CPX_SV_REF_SEGMENTS).size() > 1)
                        .cache();

        final MultiSegmentsCpxVariantExtractor multiSegmentsCpxVariantExtractor = new MultiSegmentsCpxVariantExtractor();
        final List<VariantContext> sourceWithLessAnnotations = multiSegmentCalls
                .flatMap(vc -> multiSegmentsCpxVariantExtractor.extract(vc, referenceBroadcast.getValue()).iterator()).collect();

        final List<VariantContext> sourceWithMoreAnnotations =
                reInterpretMultiSegmentComplexVarThroughAlignmentPairIteration(multiSegmentCalls,
                        svDiscoveryInputMetaData, assemblyRawAlignments);

        final List<VariantContext> reInterpretMultiSegmentsCalls = removeDuplicates(sourceWithLessAnnotations, sourceWithMoreAnnotations);
        multiSegmentCalls.unpersist(false);

        return new ExtractedSimpleVariants(reInterpretedZeroAndOneSegmentCalls, reInterpretMultiSegmentsCalls);
    }

    //==================================================================================================================

    @VisibleForTesting
    static final class RelevantAttributes implements Serializable {
        private static final long serialVersionUID = 1L;

        private final String id;
        private final List<SimpleInterval> referenceSegments;
        private final List<String> altArrangements;

        @VisibleForTesting
        RelevantAttributes(final VariantContext multiSegmentComplexVar) {
            id = multiSegmentComplexVar.getID();
            referenceSegments = SVUtils.getAttributeAsStringList(multiSegmentComplexVar, CPX_SV_REF_SEGMENTS)
                    .stream().map(SimpleInterval::new).collect(Collectors.toList());
            altArrangements = SVUtils.getAttributeAsStringList(multiSegmentComplexVar, CPX_EVENT_ALT_ARRANGEMENTS);
        }
    }

    /**
     * Send relevant contigs for re-interpretation via the pair-iteration way of scanning the alignments for interpretation.
     *
     * Re-interpret CPX vcf records whose
     * {@link org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants#CPX_SV_REF_SEGMENTS}
     * has more than one entries, aka "multi-segment" calls.
     *
     * Exist basically to extract insertions, because
     * deletions and inversions are relatively easy to be extracted by
     * {@link org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SegmentedCpxVariantSimpleVariantExtractor.MultiSegmentsCpxVariantExtractor}
     *
     * @return the {@link SimpleSVType}-d variants extracted from the input that are consistent with the description in the complex variants
     */
    private static List<VariantContext> reInterpretMultiSegmentComplexVarThroughAlignmentPairIteration(
            final JavaRDD<VariantContext> multiSegmentCalls,
            final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
            final JavaRDD<GATKRead> assemblyRawAlignments) {

        final Map<String, RelevantAttributes> contigNameToCpxVariantAttributes =
                multiSegmentCalls
                        .flatMapToPair(complex -> {
                            final RelevantAttributes relevantAttributes = new RelevantAttributes(complex);
                            return SVUtils.getAttributeAsStringList(complex, CONTIG_NAMES).stream()
                                    .map(name -> new Tuple2<>(name, relevantAttributes))
                                    .iterator();
                        })
                        .collectAsMap();

        // resend the relevant contigs through the pair-iteration-ed path
        final Set<String> relevantContigs = new HashSet<>( contigNameToCpxVariantAttributes.keySet() );
        final JavaRDD<GATKRead> relevantAlignments = assemblyRawAlignments.filter(read -> relevantContigs.contains(read.getName()));
        final JavaRDD<AlignedContig> analysisReadyContigs =
                SvDiscoverFromLocalAssemblyContigAlignmentsSpark
                        .preprocess(svDiscoveryInputMetaData, relevantAlignments)
                        .getContigsWithSignatureClassifiedAsComplex()
                        .map(AssemblyContigWithFineTunedAlignments::getSourceContig);

        List<VariantContext> pairIterationReInterpreted = ContigChimericAlignmentIterativeInterpreter
                .discoverVariantsFromChimeras(svDiscoveryInputMetaData, analysisReadyContigs);

        final Broadcast<ReferenceMultiSparkSource> referenceBroadcast = svDiscoveryInputMetaData.getReferenceData().getReferenceBroadcast();
        return pairIterationReInterpreted.stream()
                .map(vc -> {
                    final List<String> consistentComplexVariantIDs =
                            SVUtils.getAttributeAsStringList(vc, CONTIG_NAMES).stream()
                                    .map(contigNameToCpxVariantAttributes::get)
                                    .filter(attributes -> isConsistentWithCPX(vc, attributes))
                                    .map(attributes -> attributes.id)
                                    .collect(Collectors.toList());
                    if ( consistentComplexVariantIDs.isEmpty()) {
                        return null;
                    } else {
                        return new VariantContextBuilder(vc)
                                .id(vc.getID() + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                                        + CPX_DERIVED_POSTFIX_STRING)
                                .attribute(CPX_EVENT_KEY,
                                        String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, consistentComplexVariantIDs))
                                .make();
                    }
                })
                .filter(Objects::nonNull)
                .map(SegmentedCpxVariantSimpleVariantExtractor::postProcessConvertShortDupToIns)
                .flatMap(simple -> postProcessConvertReplacementToFatInsOrInsAndDel(simple, referenceBroadcast.getValue()))
                .collect(Collectors.toList());
    }

    /**
     * Convert short, i.e. duplicated range is < 50 bp, duplication call to insertion call.
     */
    @VisibleForTesting
    static VariantContext postProcessConvertShortDupToIns(final VariantContext simple) {
        final String type = simple.getAttributeAsString(SVTYPE, "");
        if ( type.equals(SimpleSVType.SupportedType.DUP.name()) ) {
            final SimpleInterval duplicatedRegion = new SimpleInterval(simple.getAttributeAsString(DUP_REPEAT_UNIT_REF_SPAN, ""));
            if (duplicatedRegion.size() > EVENT_SIZE_THRESHOLD) {
                return simple;
            } else {
                return new VariantContextBuilder(simple)
                        .alleles(Arrays.asList(simple.getReference(), altSymbAlleleIns))
                        .rmAttribute(SVTYPE)
                        .attribute(SVTYPE, SimpleSVType.SupportedType.INS.name())
                        .make();
            }
        } else
            return simple;
    }

    /**
     * Pair-iteration way of extracting simple variants reports replacement calls as a single DEL with
     * inserted sequence annotations.
     * This utility breaks that into:
     *      when the inserted sequence is long enough, an extra insertion call
     *      when the deleted range is not long enough, replace with fat insertion.
     */
    @VisibleForTesting
    static Stream<VariantContext> postProcessConvertReplacementToFatInsOrInsAndDel(final VariantContext simple,
                                                                                   final ReferenceMultiSparkSource reference) {
        final String type = simple.getAttributeAsString(SVTYPE, "");
        if ( type.equals(SimpleSVType.SupportedType.DEL.name()) ) {
            final int deletionLen = - simple.getAttributeAsInt(SVLEN, 0);
            final int insLen = simple.getAttributeAsInt(INSERTED_SEQUENCE_LENGTH, 0);
            if (insLen > EVENT_SIZE_THRESHOLD && deletionLen > EVENT_SIZE_THRESHOLD) { // case 1: insertion and deletion, linked

                final Map<String, Object> attributes = new HashMap<>( simple.getAttributes() );
                attributes.remove(INSERTED_SEQUENCE_MAPPINGS);
                attributes.remove(SVLEN);
                attributes.remove(SVTYPE);

                VariantContextBuilder newInsertion = makeInsertion(simple.getContig(), simple.getStart(), simple.getStart(), insLen, simple.getReference());
                attributes.forEach(newInsertion::attribute);
                newInsertion.rmAttribute(HOMOLOGY).rmAttribute(HOMOLOGY_LENGTH);
                newInsertion.rmAttribute(VCFConstants.END_KEY).attribute(VCFConstants.END_KEY, simple.getStart());

                VariantContextBuilder newDeletion = makeDeletion(new SimpleInterval(simple.getContig(), simple.getStart(), simple.getEnd()), simple.getReference());
                attributes.forEach(newDeletion::attribute);
                newDeletion.rmAttribute(INSERTED_SEQUENCE).rmAttribute(INSERTED_SEQUENCE_LENGTH).rmAttribute(SEQ_ALT_HAPLOTYPE);

                // cross linking
                newInsertion.attribute(LINK, makeID(SimpleSVType.SupportedType.DEL.name(), simple.getContig(), simple.getStart(), simple.getEnd()));
                newDeletion.attribute(LINK, makeID(SimpleSVType.SupportedType.INS.name(), simple.getContig(), simple.getStart(), simple.getStart()));

                return Stream.of(newDeletion.make(), newInsertion.make());
            } else if (insLen > EVENT_SIZE_THRESHOLD && deletionLen <= EVENT_SIZE_THRESHOLD) { // case 2: insertion with micro deletion
                String fatInsertionID = simple.getID().replace("DEL", "INS");
                final Map<String, Object> attributes = new HashMap<>( simple.getAttributes() );
                attributes.remove(INSERTED_SEQUENCE_MAPPINGS);
                attributes.remove(HOMOLOGY_LENGTH);
                attributes.remove(HOMOLOGY);
                attributes.remove(SVLEN);
                attributes.remove(SVTYPE);
                byte[] referenceBases = getReferenceBases(new SimpleInterval(simple.getContig(), simple.getStart(), simple.getEnd()), reference);
                VariantContextBuilder fatInsertion = makeInsertion(simple.getContig(), simple.getStart(), simple.getEnd(), insLen,
                                                                    Allele.create(referenceBases, true));
                attributes.forEach(fatInsertion::attribute);
                fatInsertion.id(fatInsertionID);
                return Stream.of(fatInsertion.make());
            } else if (insLen <= EVENT_SIZE_THRESHOLD && deletionLen > EVENT_SIZE_THRESHOLD) { // case 3:deletion with micro insertion
                return Stream.of(simple);
            } else { // case 4: neither is large enough, rare but possible
                return Stream.empty();
            }
        } else
            return Stream.of(simple);
    }

    // TODO: 3/26/18 here we check consistency only for DEL calls, and reject all INV calls (they will be extracted via MultiSegmentsCpxVariantExtractor), and INS consistency check is difficult
    /**
     * @param simple        simple variant derived from pair-iteration logic that is to be checked
     * @param attributes    source CPX variant attributes
     */
    @VisibleForTesting
    static boolean isConsistentWithCPX(final VariantContext simple,
                                       final RelevantAttributes attributes) {

        final String typeString = simple.getAttributeAsString(SVTYPE, "");

        if ( typeString.equals(SimpleSVType.SupportedType.DEL.name()) ) {
            final List<SimpleInterval> refSegments = attributes.referenceSegments;
            final List<String> altArrangement = attributes.altArrangements;

            final Tuple3<Set<SimpleInterval>, Set<Integer>, List<Integer>> missingAndPresentAndInvertedSegments =
                    getMissingAndPresentAndInvertedSegments(refSegments, altArrangement);
            final Set<SimpleInterval> missingSegments = missingAndPresentAndInvertedSegments._1();

            return deletionConsistencyCheck(simple, missingSegments);
        } else if ( typeString.equals(SimpleSVType.SupportedType.INV.name()) ) {
            return false;
        } else
            return true;
    }

    @VisibleForTesting
    static boolean deletionConsistencyCheck(final VariantContext simple, final Set<SimpleInterval> missingSegments) {
        if (missingSegments.isEmpty()) return false;

        final SimpleInterval deletedRange = new SimpleInterval(simple.getContig(), simple.getStart() + 1, simple.getEnd());
        // dummy number for chr to be used in constructing SVInterval, since 2 input AI's both map to the same chr by this point
        final int dummyChr = 0;
        final SVInterval intervalOne = new SVInterval(dummyChr, deletedRange.getStart() - 1, deletedRange.getEnd());

        for (final SimpleInterval missing : missingSegments) {
            if ( ! missing.overlaps(deletedRange) )
                return false;
            final SVInterval intervalTwo = new SVInterval(dummyChr, missing.getStart() - 1, missing.getEnd());
            // allow 1-base fuzziness from either end
            if ( Math.abs(missing.size() - deletedRange.size()) > 2 )
                return false;
            if( 2 >= Math.abs( Math.min(missing.size(), deletedRange.size()) - intervalTwo.overlapLen(intervalOne) ) ){
                return true;
            }
        }
        return false;
    }

    /**
     * Exist for equals() and hashCode()
     */
    private static final class AnnotatedInterval {

        private final VariantContext sourceVC; // NOTE: omitted in equals() and hashCode() on purpose

        final SimpleInterval interval;
        final String id;
        final String type;
        final int svlen;
        final List<Allele> alleles;

        private AnnotatedInterval(final VariantContext vc) {
            sourceVC = vc;
            interval = new SimpleInterval( vc.getContig(), vc.getStart(), vc.getEnd());
            id = vc.getID();
            type = vc.getAttributeAsString(SVTYPE, "");
            svlen = vc.getAttributeAsInt(SVLEN, 0);
            alleles = vc.getAlleles();
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final AnnotatedInterval interval1 = (AnnotatedInterval) o;

            if (svlen != interval1.svlen) return false;
            if (!interval.equals(interval1.interval)) return false;
            if (!id.equals(interval1.id)) return false;
            if (!type.equals(interval1.type)) return false;
            return alleles.equals(interval1.alleles);
        }

        @Override
        public int hashCode() {
            int result = interval.hashCode();
            result = 31 * result + id.hashCode();
            result = 31 * result + type.hashCode();
            result = 31 * result + svlen;
            result = 31 * result + alleles.hashCode();
            return result;
        }
    }

    /**
     * For constructing a map from {@link AnnotatedInterval} to source complex variant IDs and
     * their associated assembly contig names.
     */
    private static Map<AnnotatedInterval, Tuple2<TreeSet<String>, TreeSet<String>>>
    getAnnotatedIntervalToSourceCpxIDsAndContigNames(final List<VariantContext> extractedSimpleVariants) {
        // TODO: 5/11/18 this is suboptimal:
        // a round trip to AnnotatedInterval because some CPX variants themselves are duplicated,
        // i.e. their alt seq, extracted from different assembly contigs, only differ slightly.
        return extractedSimpleVariants.stream().map(AnnotatedInterval::new).collect(Collectors.toCollection(HashSet::new))
                .stream().map(ai -> ai.sourceVC)
                .collect(Collectors.toMap(AnnotatedInterval::new,
                        simpleVC -> {
                    final TreeSet<String> complexEvents = new TreeSet<>(SVUtils.getAttributeAsStringList(simpleVC, CPX_EVENT_KEY));
                    final TreeSet<String> sourceCtgNames = new TreeSet<>(SVUtils.getAttributeAsStringList(simpleVC, CONTIG_NAMES));
                    return new Tuple2<>(complexEvents, sourceCtgNames);
                })
                ); // hashMap is good enough for us
    }

    /**
     * Exist because the two ways to re-interpret simple variants via
     * {@link MultiSegmentsCpxVariantExtractor}
     * and via
     * {@link #reInterpretMultiSegmentComplexVarThroughAlignmentPairIteration(JavaRDD, SvDiscoveryInputMetaData, JavaRDD)}
     * could give essentially the same variants.
     */
    @VisibleForTesting
    static List<VariantContext> removeDuplicates(final List<VariantContext> sourceWithLessAnnotations,
                                                 final List<VariantContext> sourceWithMoreAnnotations) {

        final Map<AnnotatedInterval, Tuple2<TreeSet<String>, TreeSet<String>>> rangeToAnnotationsFromSourceWithLessAnnotations =
                getAnnotatedIntervalToSourceCpxIDsAndContigNames(sourceWithLessAnnotations);
        final Map<AnnotatedInterval, Tuple2<TreeSet<String>, TreeSet<String>>> rangeToAnnotationsFromSourceWithMoreAnnotations =
                getAnnotatedIntervalToSourceCpxIDsAndContigNames(sourceWithMoreAnnotations);

        final List<VariantContext> result = new ArrayList<>(sourceWithMoreAnnotations.size() + sourceWithLessAnnotations.size());
        for (final Map.Entry<AnnotatedInterval, Tuple2<TreeSet<String>, TreeSet<String>>> entry: rangeToAnnotationsFromSourceWithMoreAnnotations.entrySet()) {
            final AnnotatedInterval interval = entry.getKey();
            final Tuple2<TreeSet<String>, TreeSet<String>> sourceAttributes = entry.getValue();
            final Tuple2<TreeSet<String>, TreeSet<String>> anotherSourceAttributes = rangeToAnnotationsFromSourceWithLessAnnotations.get(interval);
            if (anotherSourceAttributes == null) { // variant unique to one source
                result.add( interval.sourceVC );
            } else { // found duplicate, merge annotations
                final TreeSet<String> sourceCpxIDs = sourceAttributes._1;
                final TreeSet<String> sourceCtgNames = sourceAttributes._2;
                sourceCpxIDs.addAll(anotherSourceAttributes._1);
                sourceCtgNames.addAll(anotherSourceAttributes._2);
                final VariantContextBuilder variant = new VariantContextBuilder(interval.sourceVC)
                        .rmAttribute(CPX_EVENT_KEY)
                        .attribute(CPX_EVENT_KEY, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, sourceCpxIDs))
                        .rmAttribute(CONTIG_NAMES)
                        .attribute(CONTIG_NAMES, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, sourceCtgNames));
                result.add( variant.make());

                rangeToAnnotationsFromSourceWithLessAnnotations.remove(interval); // remove from the other source
            }
        }

        // now anotherSource has only unique records
        rangeToAnnotationsFromSourceWithLessAnnotations.keySet().forEach(interval -> result.add(interval.sourceVC));

        return result;
    }

    //==================================================================================================================

    abstract List<VariantContext> extract(final VariantContext complexVC, final ReferenceMultiSparkSource reference);

    @VisibleForTesting
    static final class ZeroAndOneSegmentCpxVariantExtractor extends SegmentedCpxVariantSimpleVariantExtractor {
        private static final long serialVersionUID = 1L;

        /**
         * Depending on how the ref segment is present in alt arrangement (if at all), logic as follows (order is important):
         * <ul>
         *     <li>
         *         if ref segment appear inverted and large enough
         *         <ul>
         *             <li> INV call is warranted </li>
         *             <li> INS call(s) before and after the INV, if inserted sequence long enough </li>
         *         </ul>
         *     </li>
         *
         *     <li>
         *         otherwise if ref segment is present as-is, i.e. no deletion call can be made,
         *         make insertion calls when possible
         *     </li>
         *     <li>
         *         otherwise
         *         <ul>
         *             <li> if the segment is large enough, make a DEL call, and insertion calls when possible </li>
         *             <li> otherwise a single fat INS call</li>
         *         </ul>
         *     </li>
         * </ul>
         *
         * <p>
         *     Note that the above logic has a bias towards getting INV calls, because
         *     when the (large enough) reference segment appears both as-is and inverted,
         *     the above logic will emit at least an INV call,
         *     whereas the (inverted) duplication(s) could also be reported as an DUP call as well, but...
         * </p>
         */
        @Override
        List<VariantContext> extract(final VariantContext complexVC, final ReferenceMultiSparkSource reference) {

            final List<String> segments = SVUtils.getAttributeAsStringList(complexVC, CPX_SV_REF_SEGMENTS);
            if (segments.isEmpty()) return whenZeroSegments(complexVC, reference);

            final SimpleInterval refSegment = new SimpleInterval(segments.get(0));
            final List<String> altArrangement = SVUtils.getAttributeAsStringList(complexVC, CPX_EVENT_ALT_ARRANGEMENTS);
            final int altSeqLength = complexVC.getAttributeAsString(SEQ_ALT_HAPLOTYPE, "").length();

            final List<VariantContextBuilder> result = new ArrayList<>();

            final int asIsAppearanceIdx = altArrangement.indexOf("1");
            final int invertedAppearanceIdx = altArrangement.indexOf("-1");
            if (invertedAppearanceIdx != -1 && refSegment.size() > EVENT_SIZE_THRESHOLD) { // inversion call
                whenInversionIsWarranted(refSegment, invertedAppearanceIdx, altArrangement, reference, result);
            } else if (asIsAppearanceIdx != -1) { // no inverted appearance or appear inverted but not large enough, and in the mean time appear as-is, so no deletion
                whenNoDeletionIsAllowed(refSegment, asIsAppearanceIdx, altArrangement, altSeqLength, reference, result);
            } else { // no as-is appearance && (inverted appearance might present not not large enough)
                whenNoInvAndNoAsIsAppearance(refSegment, altSeqLength, reference, result);
            }

            final String sourceID = complexVC.getID();
            final List<String> evidenceContigs = SVUtils.getAttributeAsStringList(complexVC, CONTIG_NAMES);
            final List<String> mappingQualities = SVUtils.getAttributeAsStringList(complexVC, MAPPING_QUALITIES);
            final int maxAlignLength = complexVC.getAttributeAsInt(MAX_ALIGN_LENGTH, 0);
            return result.stream()
                    .map(vc -> vc.attribute(CPX_EVENT_KEY, sourceID).attribute(CONTIG_NAMES, evidenceContigs)
                                 .attribute(MAPPING_QUALITIES, mappingQualities)
                                 .attribute(MAX_ALIGN_LENGTH, maxAlignLength).make())
                    .collect(Collectors.toList());
        }

        private List<VariantContext> whenZeroSegments(final VariantContext complexVC, final ReferenceMultiSparkSource reference) {
            final Allele anchorBaseRefAllele = getAnchorBaseRefAllele(complexVC.getContig(), complexVC.getStart(), reference);
            final int altSeqLength = complexVC.getAttributeAsString(SEQ_ALT_HAPLOTYPE, "").length() - 2;
            final List<String> mappingQualities = SVUtils.getAttributeAsStringList(complexVC, MAPPING_QUALITIES);
            final int maxAlignLength = complexVC.getAttributeAsInt(MAX_ALIGN_LENGTH, 0);
            final VariantContext insertion = makeInsertion(complexVC.getContig(), complexVC.getStart(), complexVC.getStart(), altSeqLength, anchorBaseRefAllele)
                    .attribute(CPX_EVENT_KEY, complexVC.getID())
                    .attribute(CONTIG_NAMES, complexVC.getAttribute(CONTIG_NAMES))
                    .attribute(MAPPING_QUALITIES, mappingQualities)
                    .attribute(MAX_ALIGN_LENGTH, maxAlignLength)
                    .make();
            return Collections.singletonList(insertion);
        }

        private static void whenInversionIsWarranted(final SimpleInterval refSegment, final int invertedAppearanceIdx,
                                                     final List<String> altArrangement, final ReferenceMultiSparkSource reference,
                                                     final List<VariantContextBuilder> result) {

            final Allele anchorBaseRefAllele = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getStart(), reference);
            result.add( makeInversion(refSegment, anchorBaseRefAllele) );

            // further check if alt seq length is long enough to trigger an insertion as well,
            // but guard against case smallIns1 + INV + smallIns2, in theory one could annotate the inversion
            // with micro-insertions if that's the case, but we try to have minimal annotations here
            final Allele anchorBaseRefAlleleFront = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getStart() - 1, reference);
            final Allele anchorBaseRefAlleleRear  = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getEnd(), reference);
            extractFrontAndRearInsertions(refSegment, invertedAppearanceIdx, altArrangement,
                    anchorBaseRefAlleleFront, anchorBaseRefAlleleRear, result);
        }

        private static void whenNoDeletionIsAllowed( final SimpleInterval refSegment, final int asIsAppearanceIdx,
                                                     final List<String> altArrangement, final int altSeqLength,
                                                     final ReferenceMultiSparkSource reference, final List<VariantContextBuilder> result) {
            final int segmentSize = refSegment.size();
            if (altSeqLength - segmentSize > EVENT_SIZE_THRESHOLD ) { // long enough net gain to trigger insertion calls
                // distinguish between cases {"1", ....}, {....., "1"}, and {....., "1", ....} to know where to place the insertion
                final Allele anchorBaseRefAlleleFront = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getStart() - 1, reference);
                final Allele anchorBaseRefAlleleRear = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getEnd(), reference);
                if ( altArrangement.get(altArrangement.size() - 1).equals("1") ) { // {....., "1"} -> front insertion
                    final VariantContextBuilder frontIns =
                            SegmentedCpxVariantSimpleVariantExtractor.makeInsertion(refSegment.getContig(),
                                    refSegment.getStart() - 1, refSegment.getStart() - 1,
                                    altSeqLength - segmentSize, anchorBaseRefAlleleFront);
                    result.add(frontIns);
                } else if ( altArrangement.get(0).equals("1") ) { // {"1", ....} -> rear insertion
                    final VariantContextBuilder rearIns =
                            SegmentedCpxVariantSimpleVariantExtractor.makeInsertion(refSegment.getContig(),
                                    refSegment.getEnd(), refSegment.getEnd(),
                                    altSeqLength - segmentSize, anchorBaseRefAlleleFront);
                    result.add(rearIns);
                } else { // {....., "1", ....} -> collect new insertion length before and after
                    extractFrontAndRearInsertions(refSegment, asIsAppearanceIdx, altArrangement,
                            anchorBaseRefAlleleFront, anchorBaseRefAlleleRear, result);
                }
            }
        }

        private static void whenNoInvAndNoAsIsAppearance( final SimpleInterval refSegment, final int altSeqLength,
                                                          final ReferenceMultiSparkSource reference, final List<VariantContextBuilder> result) {
            if ( refSegment.size() > EVENT_SIZE_THRESHOLD ) { // a deletion call must be present

                final Allele anchorBaseRefAlleleFront = getAnchorBaseRefAllele(refSegment.getContig(), refSegment.getStart(), reference);

                // need left shift because the segment boundaries are shared by REF and ALT
                result.add( makeDeletion(new SimpleInterval(refSegment.getContig(), refSegment.getStart(), refSegment.getEnd() - 1),
                                         anchorBaseRefAlleleFront) );

                // if the replacing sequence is long enough to trigger an insertion as well
                if (altSeqLength - 2 > EVENT_SIZE_THRESHOLD) {
                    result.add(makeInsertion(refSegment.getContig(), refSegment.getStart(), refSegment.getStart(), altSeqLength, anchorBaseRefAlleleFront));
                }
            } else if ( altSeqLength - 2 > EVENT_SIZE_THRESHOLD ){ // ref segment not long enough to merit an INV or DEL, so a fat INS, if size is enough
                final Allele fatInsertionRefAllele =
                        Allele.create(getReferenceBases(new SimpleInterval(refSegment.getContig(), refSegment.getStart(), refSegment.getEnd() - 1), reference), true);
                result.add( makeInsertion(refSegment.getContig(), refSegment.getStart(), refSegment.getEnd() - 1,
                        altSeqLength - refSegment.size(), fatInsertionRefAllele) );
            }
        }

        private static void extractFrontAndRearInsertions(final SimpleInterval refSegment, final int segmentIdx,
                                                          final List<String> altArrangement,
                                                          final Allele anchorBaseRefAlleleFront,
                                                          final Allele anchorBaseRefAlleleRear,
                                                          final List<VariantContextBuilder> result) {

            final List<Integer> segmentLen = Collections.singletonList(refSegment.size());

            final SimpleInterval frontInsPos = SVUtils.makeOneBpInterval(refSegment.getContig(), refSegment.getStart() - 1);
            final VariantContextBuilder frontIns =
                    getInsFromOneEnd(true, segmentIdx, frontInsPos, anchorBaseRefAlleleFront, segmentLen, altArrangement, true);
            if (frontIns != null)
                result.add(frontIns);

            final SimpleInterval rearInsPos = SVUtils.makeOneBpInterval(refSegment.getContig(), refSegment.getEnd());
            final VariantContextBuilder rearIns =
                    getInsFromOneEnd(false, segmentIdx, rearInsPos, anchorBaseRefAlleleRear, segmentLen, altArrangement, true);
            if (rearIns != null)
                result.add(rearIns);
        }
    }

    @VisibleForTesting
    static final class MultiSegmentsCpxVariantExtractor extends SegmentedCpxVariantSimpleVariantExtractor {
        private static final long serialVersionUID = 1L;

        @Override
        List<VariantContext> extract(final VariantContext complexVC, final ReferenceMultiSparkSource reference) {

            final List<SimpleInterval> refSegments =
                    SVUtils.getAttributeAsStringList(complexVC, CPX_SV_REF_SEGMENTS).stream()
                            .map(SimpleInterval::new)
                            .collect(Collectors.toList());

            final List<String> altArrangement = SVUtils.getAttributeAsStringList(complexVC, CPX_EVENT_ALT_ARRANGEMENTS);

            final Tuple3<Set<SimpleInterval>, Set<Integer>, List<Integer>> missingAndPresentAndInvertedSegments = getMissingAndPresentAndInvertedSegments(refSegments, altArrangement);
            final Set<SimpleInterval> missingSegments = missingAndPresentAndInvertedSegments._1();
            final Set<Integer> presentSegments = missingAndPresentAndInvertedSegments._2();
            final List<Integer> invertedSegments = missingAndPresentAndInvertedSegments._3();

            final List<VariantContextBuilder> result = new ArrayList<>();

            // if affected ref sequence found as is (trusting the aligner), then only output front and/or rear insertions
            final int idx = findAllSegments(altArrangement, refSegments.size());
            if ( idx >= 0 ) {
                whenAllSegmentsAppearAsIs(complexVC, reference, refSegments, altArrangement, result, idx);
            } else {

                // inversions
                if (!invertedSegments.isEmpty()) {
                    extractInversions(reference, refSegments, presentSegments, invertedSegments, result);
                }

                // deletions
                if (!missingSegments.isEmpty()) {
                    extractDeletions(reference, missingSegments, result);
                }

                // head and tail insertions only
                extractFrontAndRearInsertions(complexVC, refSegments, altArrangement, reference, result);
            }

            final String sourceID = complexVC.getID();
            final List<String> evidenceContigs = SVUtils.getAttributeAsStringList(complexVC, CONTIG_NAMES);
            final List<String> mappingQualities = SVUtils.getAttributeAsStringList(complexVC, MAPPING_QUALITIES);
            final int maxAlignLength = complexVC.getAttributeAsInt(MAX_ALIGN_LENGTH, 0);

            return result.stream()
                    .map(vc -> vc.attribute(CPX_EVENT_KEY, sourceID).attribute(CONTIG_NAMES, evidenceContigs)
                                 .attribute(MAPPING_QUALITIES, mappingQualities)
                                 .attribute(MAX_ALIGN_LENGTH, maxAlignLength).make())
                    .collect(Collectors.toList());
        }

        /**
         * Given {@code altArrangement} and count of segments, return the index in {@code altArrangement}
         * pointing to "1" where all segments contiguously appear after that, i.e. the affected reference region
         * appear as is (with other insertions, duplications, etc at either end) according to {@code altArrangement}
         *
         * Example:
         *  ......, 1, 2, 3, 4, .....
         *  with 4 segments,
         *  and the index of the "1" is 2,
         *  then this function returns 2.
         *  but if the altArrangement is
         *  ......, 1, 2, 3, , .....
         *  the function returns -1 because not all segments appear as-is in the description.
         */
        @VisibleForTesting
        static int findAllSegments(final List<String> altArrangement, final int segmentCount) {
            int idx = -1;
            int currentlyLookingForSegment = segmentCount;
            final String segmentCountString = String.valueOf(segmentCount);
            for (int i = altArrangement.size() - 1; i >= 0 ; --i) { // reversely because we want to follow left-justify convention
                final String description = altArrangement.get(i);
                if ( description.equals( String.valueOf(currentlyLookingForSegment) ) ) {
                    if (currentlyLookingForSegment == 1) return i;
                    --currentlyLookingForSegment;
                } else {
                    currentlyLookingForSegment = description.equals( segmentCountString ) ? segmentCount - 1 : segmentCount;
                    idx = -1;
                }
            }
            return idx;
        }

        private static void whenAllSegmentsAppearAsIs(final VariantContext complexVC, final ReferenceMultiSparkSource reference,
                                                      final List<SimpleInterval> refSegments, final List<String> altArrangement,
                                                      final List<VariantContextBuilder> result, final int idx) {
            final List<Integer> refSegmentLengths = refSegments.stream().map(SimpleInterval::size).collect(Collectors.toList());

            if ( idx != 0 ) { // e.g. 4 segments, and alt arrangement is ......, 1,2,3,4, there could be (that is, if long enough) front insertion
                final SimpleInterval insertionPos = new SimpleInterval(complexVC.getContig(),
                        complexVC.getStart() - 1, complexVC.getStart() - 1);
                final Allele anchorBaseRefAlleleFront = getAnchorBaseRefAllele(insertionPos.getContig(), insertionPos.getStart(), reference);
                final VariantContextBuilder frontIns = getInsFromOneEnd(true, idx, insertionPos,
                        anchorBaseRefAlleleFront, refSegmentLengths, altArrangement, true);
                if (frontIns != null) result.add(frontIns);
            }
            if ( idx + refSegments.size() - 1 < altArrangement.size() - 1 ) { // e.g. there's more after 1,2,3,4,..., there could be (that is, if long enough) front insertion
                final SimpleInterval insertionPos = new SimpleInterval(complexVC.getContig(), complexVC.getEnd(), complexVC.getEnd());
                final byte[] refBases = getReferenceBases(insertionPos, reference);
                final Allele anchorBaseRefAlleleRear = Allele.create(refBases, true);
                final VariantContextBuilder rearIns = getInsFromOneEnd(false, idx + refSegments.size() - 1, insertionPos,
                        anchorBaseRefAlleleRear, refSegmentLengths, altArrangement, true);
                if (rearIns != null) result.add(rearIns);
            }
        }

        private void extractInversions( final ReferenceMultiSparkSource reference, final List<SimpleInterval> refSegmentIntervals,
                                        final Set<Integer> presentSegments, final List<Integer> invertedSegments,
                                        final List<VariantContextBuilder> result) {
            final List<VariantContextBuilder> inversions =
                    invertedSegments.stream()
                            // large enough; in addition, if both as-is and inverted versions exist, treat as insertions instead of inversions: unlike 1-segment calls, where we don't have consistency problems
                        .filter(i -> refSegmentIntervals.get(i - 1).size() > EVENT_SIZE_THRESHOLD && (!presentSegments.contains(i)))
                        .map(i -> {
                            final SimpleInterval invertedSegment = refSegmentIntervals.get(i - 1);
                            final byte[] ref = getReferenceBases(SVUtils.makeOneBpInterval(invertedSegment.getContig(), invertedSegment.getStart()), reference);
                            final Allele refAllele = Allele.create(ref, true);
                            return makeInversion(invertedSegment, refAllele);
                        })
                        .collect(Collectors.toList());
            result.addAll(inversions);
        }

        private void extractDeletions( final ReferenceMultiSparkSource reference, final Set<SimpleInterval> missingSegments,
                                       final List<VariantContextBuilder> result) {
            final List<VariantContextBuilder> deletions = compactifyMissingSegments(missingSegments).stream()
                    .filter(gone -> gone.size() > EVENT_SIZE_THRESHOLD) // large enough
                    .map(gone -> {
                        final byte[] ref = getReferenceBases(SVUtils.makeOneBpInterval(gone.getContig(), gone.getStart()), reference);
                        final Allele refAllele = Allele.create(ref, true);
                        return makeDeletion(new SimpleInterval(gone.getContig(), gone.getStart(), gone.getEnd() - 1), refAllele);
                    })
                    .collect(Collectors.toList());
            result.addAll(deletions);
        }
        /**
         * Compactify missingSegments for case when two neighboring segments are both gone, to avoid cases when
         * 1) neither segment is large enough
         * 2) calling two small deletions while one should call a big deletion
         */
        @VisibleForTesting
        static List<SimpleInterval> compactifyMissingSegments(final Set<SimpleInterval> missingSegments) {
            if (missingSegments.size() == 1)
                return Collections.singletonList(missingSegments.iterator().next());

            // first sort
            final List<SimpleInterval> sortedMissingSegments = missingSegments.stream()
                    .sorted(Comparator.comparing(SimpleInterval::getStart)) // two segments will NEVER have the same start or overlap on more than one base
                    .collect(Collectors.toList());
            final List<SimpleInterval> result = new ArrayList<>(missingSegments.size());
            Iterator<SimpleInterval> iterator = sortedMissingSegments.iterator();
            SimpleInterval current = iterator.next();
            while (iterator.hasNext()) {
                SimpleInterval next = iterator.next();
                if (current.overlapsWithMargin(next, 1)) {
                    current = new SimpleInterval(current.getContig(), current.getStart(), next.getEnd());
                } else {
                    result.add(current);
                    current = next;
                }
            }
            result.add(current);
            return result;
        }

        private void extractFrontAndRearInsertions(final VariantContext complexVC, final List<SimpleInterval> refSegmentIntervals,
                                                   final List<String> altArrangement, final ReferenceMultiSparkSource reference,
                                                   final List<VariantContextBuilder> result) {
            final List<Integer> refSegmentLengths = refSegmentIntervals.stream().map(SimpleInterval::size).collect(Collectors.toList());
            // index pointing to first appearance of ref segment (inverted or not) in altArrangement, from either side
            int firstRefSegmentIdx = 0; // first front
            for (final String description : altArrangement) {
                if ( descriptionIndicatesInsertion(description)) {
                    ++firstRefSegmentIdx;
                } else {
                    break;
                }
            }
            if (firstRefSegmentIdx > 0) {
                final SimpleInterval startAndStop = SVUtils.makeOneBpInterval(complexVC.getContig(), complexVC.getStart());
                final Allele anchorBaseRefAlleleFront = Allele.create(getReferenceBases(startAndStop, reference), true);
                final VariantContextBuilder frontIns = getInsFromOneEnd(true, firstRefSegmentIdx, startAndStop, anchorBaseRefAlleleFront, refSegmentLengths, altArrangement, true);
                if (frontIns != null) result.add( frontIns );
            }

            firstRefSegmentIdx = altArrangement.size() - 1; // then end
            for (int i = altArrangement.size() - 1; i > -1 ; --i) {
                if ( descriptionIndicatesInsertion(altArrangement.get(i))) {
                    --firstRefSegmentIdx;
                } else {
                    break;
                }
            }

            if (firstRefSegmentIdx != altArrangement.size() - 1) {
                final int pos = complexVC.getEnd();
                final SimpleInterval insertionPos = SVUtils.makeOneBpInterval(complexVC.getContig(), pos);
                final Allele anchorBaseRefAlleleRear = Allele.create(getReferenceBases(insertionPos, reference), true);
                final VariantContextBuilder rearIns = getInsFromOneEnd(false, firstRefSegmentIdx, insertionPos, anchorBaseRefAlleleRear, refSegmentLengths, altArrangement, true);
                if (rearIns != null) result.add( rearIns );
            }
        }

        @VisibleForTesting
        static boolean descriptionIndicatesInsertion(final String description) {
            if (description.startsWith(CpxVariantCanonicalRepresentation.UNMAPPED_INSERTION))
                return true;
            return !NumberUtils.isCreatable(description); // "(-)?[0-9]+" is describing segments, we don't count them as insertions
        }
    }

    //==================================================================================================================

    /**
     * Reason for requesting increment by 1 via {@code shouldIncreaseInsLenByOne}:
     * when getting insertion length from either end,
     * there could be, but not always, a one-bp overlap between the head alignment and
     * the next alignment that continues the flow (which is not necessarily the 2nd alignment);
     * so when there's such 1-bp overlap, the insertion length should count this 1-bp overlap.
     * todo: currently all known calling code provide {@code true}, which is technically wrong, but we need alignment information for tell when to provide true/false
     *
     * @return {@code null} if the inserted sequence from the requested end is not over {@link #EVENT_SIZE_THRESHOLD}
     */
    @VisibleForTesting
    static VariantContextBuilder getInsFromOneEnd(final boolean fromFront, final int idxFirstMatch,
                                                  final SimpleInterval insertionStartAndStop, final Allele anchorBaseRefAllele,
                                                  final List<Integer> refSegmentLengths, final List<String> altArrangement,
                                                  final boolean shouldIncreaseInsLenByOne) {
        int insLen = 0;
        if (fromFront) {
            for (int i = 0; i < idxFirstMatch; ++i) {
                insLen += getInsLen( altArrangement.get(i), refSegmentLengths );
            }
        } else {
            for (int i = idxFirstMatch + 1; i < altArrangement.size(); ++i) {
                insLen += getInsLen( altArrangement.get(i), refSegmentLengths );
            }
        }

        if (shouldIncreaseInsLenByOne) ++insLen;

        if (insLen > EVENT_SIZE_THRESHOLD)
            return makeInsertion(insertionStartAndStop.getContig(), insertionStartAndStop.getStart(), insertionStartAndStop.getEnd(), insLen, anchorBaseRefAllele);
        else
            return null;
    }

    @VisibleForTesting
    static int getInsLen(final String description, final List<Integer> refSegmentLengths) {
        if (description.startsWith(CpxVariantCanonicalRepresentation.UNMAPPED_INSERTION)) {
            return Integer.valueOf(description.substring(CpxVariantCanonicalRepresentation.UNMAPPED_INSERTION.length() + 1));
        } else if ( NumberUtils.isCreatable(description) ){
            final int offset = description.startsWith("-") ? 1 : 0;
            return refSegmentLengths.get( Integer.valueOf(description.substring(offset)) - 1);
        } else {
            final int offset = description.startsWith("-") ? 1 : 0;
            return new SimpleInterval(description.substring(offset)).size();
        }
    }

    /**
     * Retrieves from the provide {@code complexVC}, reference segments described in
     * {@link org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants#CPX_SV_REF_SEGMENTS}, that are
     *   a) absent
     *   b) present as is, i.e. not inverted
     *   c) inverted
     */
    @VisibleForTesting
    static Tuple3<Set<SimpleInterval>, Set<Integer>, List<Integer>> getMissingAndPresentAndInvertedSegments(final List<SimpleInterval> refSegments,
                                                                                                            final List<String> altArrangements ) {

        final List<Integer> invertedSegments = new ArrayList<>();
        final Set<Integer> presentSegments = new TreeSet<>();
        altArrangements
                .forEach(s -> {
                    if ( s.startsWith("-") && ( !s.contains(":") )) { // some segment inverted
                        invertedSegments.add( Integer.valueOf(s.substring(1)) );
                    }
                    if ( !s.contains(":") && !s.startsWith(CpxVariantCanonicalRepresentation.UNMAPPED_INSERTION) && !s.startsWith("-") ) { // a ref segment, but not inverted
                        presentSegments.add(Integer.valueOf(s));
                    }
                });

        final Set<SimpleInterval> missingSegments = IntStream.rangeClosed(1, refSegments.size()).boxed()
                .filter(i -> !presentSegments.contains(i) && !invertedSegments.contains(i))
                .map(i -> refSegments.get(i-1))
                .collect(Collectors.toSet());

        return new Tuple3<>(missingSegments, presentSegments, invertedSegments);
    }

    // boiler-plate code block =========================================================================================

    private static Allele getAnchorBaseRefAllele(final String chr, final int pos, final ReferenceMultiSparkSource reference) {
        return Allele.create(getReferenceBases(SVUtils.makeOneBpInterval(chr, pos), reference), true);
    }

    // try not to have many try's
    static byte[] getReferenceBases(final SimpleInterval interval, final ReferenceMultiSparkSource reference) {
        try {
            return reference.getReferenceBases(interval).getBases();
        } catch (final IOException ioex) {
            throw new GATKException("Failed to extract reference bases on:" + interval, ioex);
        }
    }

    private static final Allele altSymbAlleleDel = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_DEL));
    private static final Allele altSymbAlleleIns = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_INS));
    private static final Allele altSymbAlleleInv = Allele.create(SimpleSVType.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_INV));

    /**
     * Note that {@code delRange} is expected to be pre-process to VCF spec compatible,
     * e.g. if chr1:101-200 is deleted, then {@code delRange} should be chr1:100-200
     * @param delRange
     */
    @VisibleForTesting
    static VariantContextBuilder makeDeletion(final SimpleInterval delRange, final Allele refAllele) {

        return new VariantContextBuilder()
                .chr(delRange.getContig()).start(delRange.getStart()).stop(delRange.getEnd())
                .alleles(Arrays.asList(refAllele, altSymbAlleleDel))
                .id(makeID(SimpleSVType.SupportedType.DEL.name(), delRange.getContig(), delRange.getStart(), delRange.getEnd()))
                .attribute(VCFConstants.END_KEY, delRange.getEnd())
                .attribute(SVLEN, - delRange.size() + 1)
                .attribute(SVTYPE, SimpleSVType.SupportedType.DEL.name());
    }

    @VisibleForTesting
    static VariantContextBuilder makeInsertion(final String chr, final int pos, final int end, final int svLen,
                                               final Allele refAllele) {

        return new VariantContextBuilder().chr(chr).start(pos).stop(end)
                .alleles(Arrays.asList(refAllele, altSymbAlleleIns))
                .id(makeID(SimpleSVType.SupportedType.INS.name(), chr, pos, end))
                .attribute(VCFConstants.END_KEY, end)
                .attribute(SVLEN, svLen)
                .attribute(SVTYPE, SimpleSVType.SupportedType.INS.name());
    }

    @VisibleForTesting
    static VariantContextBuilder makeInversion(final SimpleInterval invertedRegion, final Allele refAllele) {
        return new VariantContextBuilder()
                .chr(invertedRegion.getContig()).start(invertedRegion.getStart() - 1).stop(invertedRegion.getEnd())     // TODO: 5/2/18 VCF spec doesn't requst left shift by 1 for inversion POS
                .alleles(Arrays.asList(refAllele, altSymbAlleleInv))
                .id(makeID(SimpleSVType.SupportedType.INV.name(), invertedRegion.getContig(), invertedRegion.getStart() - 1, invertedRegion.getEnd()))
                .attribute(VCFConstants.END_KEY, invertedRegion.getEnd())
                .attribute(SVLEN, 0)                                                                 // TODO: 5/2/18 this is following VCF spec,
                .attribute(SVTYPE, SimpleSVType.SupportedType.INV.name());
    }
}
