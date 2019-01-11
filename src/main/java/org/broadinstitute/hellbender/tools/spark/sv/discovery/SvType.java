package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.EnumUtils;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR;


/**
 * Various types of structural variations.
 * Holding minimum information on VariantContext
 * it can hold:
 * CHR, POS, ID, REF, ALT, SVTYPE, END (optional), SVLEN(optional)
 */
public abstract class SvType {
    public static final int NO_APPLICABLE_END = -1;
    public static final int NO_APPLICABLE_LEN = -1;

    // fields whose name starts with "variant" are following the VCF spec, i.e. they can be used directly to construct VariantContextBuilder
    protected final String variantChromosome;
    protected final int variantStart;
    protected final int variantStop;
    protected final String variantId;
    protected final Allele refAllele;
    protected final Allele altAllele;
    protected final int svLen;
    protected final Map<String, Object> typeSpecificAttributes;

    protected static final Map<String, Object> noExtraAttributes = Collections.emptyMap();

    public SvType(final String variantChromosome, final int variantStart, final int variantStop, final String variantId,
                  final Allele refAllele, final Allele altAllele, final int svLen, final Map<String, Object> typeSpecificAttributes) {
        this.variantChromosome = variantChromosome;
        this.variantStart = variantStart;
        this.variantStop = variantStop;
        this.variantId = variantId;
        this.refAllele = refAllele;
        this.altAllele = altAllele;
        this.svLen = svLen;
        this.typeSpecificAttributes = typeSpecificAttributes;
    }

    public final VariantContextBuilder getBasicInformation() {
        if ( ! hasApplicableEnd() ) {
            VariantContextBuilder builder = new VariantContextBuilder()
                    .chr(variantChromosome).start(variantStart).stop(variantStart)
                    .id(variantId)
                    .alleles(new ArrayList<>(Arrays.asList(refAllele, altAllele)))
                    .attribute(GATKSVVCFConstants.SVTYPE, toString());
            typeSpecificAttributes.forEach(builder::attribute);
            return builder;
        } else { // assuming if there's valid END, there must be valid SVLEN
            VariantContextBuilder builder = new VariantContextBuilder()
                    .chr(variantChromosome).start(variantStart).stop(variantStop)
                    .id(variantId)
                    .alleles(new ArrayList<>(Arrays.asList(refAllele, altAllele)))
                    .attribute(VCFConstants.END_KEY, variantStop)
                    .attribute(GATKSVVCFConstants.SVTYPE, toString())
                    .attribute(GATKSVVCFConstants.SVLEN, svLen);
            typeSpecificAttributes.forEach(builder::attribute);
            return builder;
        }
    }

    public String getVariantChromosome() {
        return variantChromosome;
    }
    public int getVariantStart() {
        return variantStart;
    }
    public int getVariantStop() {
        return variantStop;
    }
    public final String getInternalVariantId() {
        return variantId;
    }
    Allele getRefAllele() {
        return refAllele;
    }
    public final Allele getAltAllele() {
        return altAllele;
    }
    public final int getSVLength() {
        return svLen;
    }
    public final Map<String, Object> getTypeSpecificAttributes() {
        return typeSpecificAttributes;
    }

    public abstract boolean hasApplicableEnd();
    public abstract boolean hasApplicableLength();

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SvType svType = (SvType) o;

        if (svLen != svType.svLen) return false;
        if (!variantId.equals(svType.variantId)) return false;
        if (!altAllele.equals(svType.altAllele)) return false;
        if (typeSpecificAttributes.size() != svType.typeSpecificAttributes.size()) return false;

        for ( final Map.Entry<String, Object> act : typeSpecificAttributes.entrySet() ) {
            if ( svType.typeSpecificAttributes.containsKey(act.getKey())) {
                if ( ! act.getValue().equals(svType.typeSpecificAttributes.get(act.getKey())))
                    return false;
            } else
                return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int result = variantId.hashCode();
        result = 31 * result + altAllele.hashCode();
        result = 31 * result + svLen;
        result = 31 * result + typeSpecificAttributes.hashCode();
        return result;
    }

    // TODO: 5/23/18 any better way to do this?
    public static SortedSet<String> getKnownTypes() {
        final SortedSet<String> knownTypes = new TreeSet<>( EnumUtils.getEnumMap(SimpleSVType.SupportedType.class).keySet() );

        knownTypes.add(GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR);
        knownTypes.add(GATKSVVCFConstants.BREAKEND_STR);

        return Collections.unmodifiableSortedSet(knownTypes);
    }

    public static String makeLocationString(final String chr1, final int pos1, final String chr2, final int pos2) {
        return chr1 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + pos1 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR
                + (chr2.equals(chr1) ? "" : chr2 + INTERVAL_VARIANT_ID_FIELD_SEPARATOR)
                + pos2;
    }

    public static String makeLocationString(final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype) {
        String leftContig = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getContig();
        String rightContig = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getContig();
        int pos1 = novelAdjacencyAndAltHaplotype.getLeftJustifiedLeftRefLoc().getStart();
        int pos2 = novelAdjacencyAndAltHaplotype.getLeftJustifiedRightRefLoc().getEnd();
        return makeLocationString(leftContig, pos1, rightContig, pos2);
    }

    static byte[] extractRefBases(final SimpleInterval interval, final ReferenceMultiSparkSource reference) {
        try {
            return reference.getReferenceBases(interval).getBases();
        } catch (final IOException ioex) {
            throw new GATKException("Failed to extract bases from region: " + interval.toString());
        }
    }
}
