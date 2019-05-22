package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.List;


/**
 * Converts an annotated interval representing a segment to a variant context.
 */
public class AnnotatedIntervalToSegmentVariantContextConverter {
    private AnnotatedIntervalToSegmentVariantContextConverter() {}

    /**
     * Field names that will be accepted as a call
     */
    final static List<String> callAnnotationNames = Arrays.asList(
            "CALL","Segment_Call","Call"
    );

    public final static Allele COPY_NEUTRAL_ALLELE = Allele.create(SimpleSVType.createBracketedSymbAlleleString("COPY_NEUTRAL"));
    public final static String COPY_NEUTRAL_ALLELE_STRING = COPY_NEUTRAL_ALLELE.getDisplayString();

    /**
     * Convert a segment (as annotated interval) to a {@link VariantContext}.  The variant context will have a ref allele
     *  of the starting base of the segment.  The end position will be in the attribute {@link VCFConstants#END_KEY}.
     *  Any remaining annotations in the annotated interval will be converted to string attributes.  The segment call
     *  will be the alternate allele.  Currently, the alternate allele will be one of neutral, amplification, or deletion.
     *
     * The variant context will use symbolic alleles for the calls (copy neutral, amplification, or deletion).
     * See {@link CalledCopyRatioSegment.Call} for insertions and and deletions.
     *
     * @param segment Never {@code null}
     * @param referenceContext Never {@code null}
     * @return Never {@code null}
     */
    public static VariantContext convert(final AnnotatedInterval segment, final ReferenceContext referenceContext) {
        Utils.nonNull(segment);
        Utils.nonNull(referenceContext);
        final SimpleInterval startPointInterval = new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
        final Allele refAlleleAtStart = Allele.create(referenceContext.getBases(startPointInterval), true);

        final CalledCopyRatioSegment.Call call = retrieveCall(segment);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder().chr(segment.getContig())
                .start(segment.getStart())
                .stop(segment.getEnd())
                .alleles(Arrays.asList(
                        Allele.create(refAlleleAtStart, false),
                        convertActualSegCallToAllele(call)
                ))
                .attribute(VCFConstants.END_KEY, segment.getEnd());

        // Grab every field in the annotated interval and make it a string attribute
        segment.getAnnotations().keySet().stream().forEach(k -> variantContextBuilder.attribute(k, segment.getAnnotationValue(k)));

        return variantContextBuilder
                .make();
    }

    // Returns null if no call is found in the annotated interval.
    private static CalledCopyRatioSegment.Call retrieveCall(final AnnotatedInterval segment) {
        for (final String callAnnotationName : callAnnotationNames) {
            if (segment.hasAnnotation(callAnnotationName)) {
                final CalledCopyRatioSegment.Call call = Arrays.stream(CalledCopyRatioSegment.Call.values())
                        .filter(c -> c.getOutputString().equals(segment.getAnnotationValue(callAnnotationName))).findFirst().orElse(null);
                return call;
            }
        }

        return null;
    }

    private static Allele convertActualSegCallToAllele(final CalledCopyRatioSegment.Call call) {
        if (call == null) {
            return Allele.UNSPECIFIED_ALTERNATE_ALLELE;
        }
        switch (call) {
            case DELETION:
                return Allele.create(SimpleSVType.createBracketedSymbAlleleString("DEL"), false);
            case AMPLIFICATION:
                return Allele.create(SimpleSVType.createBracketedSymbAlleleString("INS"), false);
            case NEUTRAL:
                return COPY_NEUTRAL_ALLELE;
            default:
                throw new GATKException.ShouldNeverReachHereException(call.getOutputString() + " is not represented in conversion to variant context.");
        }
    }
}
