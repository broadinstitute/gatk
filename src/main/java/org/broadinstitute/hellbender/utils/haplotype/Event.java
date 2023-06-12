package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Very simple class wrapping VariantContext when we want to be explicit that a variant is biallelic, such as
 * in a haplotype EventMap
 */
public class Event implements Locatable {
    private final String contig;
    private final int start;
    private final int stop;
    private final Allele refAllele;
    private final Allele altAllele;

    private Map<String, String> attributesForVariantContext = null;

    private static final long serialVersionUID = 1L;

    public Event(final String contig, final int start, final Allele ref, final Allele alt) {
        // TODO: Should this instead silently trim to the minimal representation?
        Utils.validateArg(ref.length() == 1 || alt.length() == 1 || differentLastBase(ref, alt), "Ref and alt alleles have same last base, hence this event is not in its minimal representation.");
        Utils.validateArg(ref.isReference(), "ref is not ref");
        this.contig = contig;
        this.start = start;
        stop = start + ref.length() - 1;
        refAllele = ref;
        altAllele = alt;
    }

    public static Event ofWithoutAttributes(final VariantContext vc) {
        Utils.validateArg(vc.isBiallelic(), "variant must be biallelic");
        return new Event(vc.getContig(), vc.getStart(), vc.getReference(), vc.getAlternateAllele(0));
    }

    // This should only be used once in the lifecycle of an event: when we make the jump from discovered event to variant context for output
    public VariantContext convertToVariantContext(final String source) {
        final VariantContext result = new VariantContextBuilder(source, contig, start, stop, Arrays.asList(refAllele, altAllele)).make();
        if (attributesForVariantContext != null) {
            attributesForVariantContext.forEach((key, value) -> result.getCommonInfo().putAttribute(key, value));
        }
        return result;
    }

    @Override
    public String getContig() { return contig; }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() { return stop; }

    public Allele refAllele() {
        return refAllele;
    }

    public Allele altAllele() {
        return altAllele;
    }

    public boolean isSNP() { return refAllele.length() == 1 && refAllele.length() == altAllele.length(); }

    public boolean isIndel() { return refAllele.length() != altAllele.length() && !altAllele.isSymbolic(); }

    public boolean isSimpleInsertion() { return refAllele.length() == 1 && altAllele.length() > 1; }

    public boolean isSimpleDeletion() { return refAllele.length() > 1 && altAllele.length() == 1; }

    public boolean isMNP() { return refAllele.length() > 1 && refAllele.length() == altAllele.length(); }

    public void setVariantAttribute(final String key, final String value) {
        if (attributesForVariantContext == null) {
            attributesForVariantContext = new HashMap<>();
        }
        attributesForVariantContext.put(key, value);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || obj.getClass() != this.getClass()) {
            return false;
        }

        final Event other = (Event) obj;

        return this.start == other.start && this.refAllele.equals(other.refAllele) && this.altAllele.equals(other.altAllele);
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder().append(start).append(refAllele).append(altAllele).hashCode();
    }

    private static boolean differentLastBase(final Allele ref, final Allele alt) {
        final byte[] refBases = ref.getBases();
        final byte[] altBases = alt.getBases();
        return refBases.length == 0 || altBases.length == 0 || refBases[refBases.length-1] != altBases[altBases.length-1];
    }

    @Override
    public String toString() {
        return "[Event @ " + contig + ":" + (start - stop == 0 ? start : start + "-" + stop) + ", " +
                refAllele + "->" + altAllele + ", " + attributesForVariantContext + "]";
    }
}
