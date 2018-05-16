package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public abstract class SVLocalContext extends VariantContext {
    public static final long serialVersionUID = 1L;

    protected SVLocalContext(final VariantContext vc) {
        super(vc);
    }

    public static boolean indicatesInversion(final VariantContext variantContext) {
        return variantContext.getAttributeAsBoolean(INV33, false)
                ||
                variantContext.getAttributeAsBoolean(INV55, false);
    }

    private static String parseAssemblyID(final String localAssemblyContigName) {
        if ( localAssemblyContigName.length() != 18
                ||
                ! localAssemblyContigName.startsWith("asm")
                ||
                ! localAssemblyContigName.contains("tig"))
            throw new IllegalArgumentException("Badly formatted contig name: " + localAssemblyContigName);

        final String assemblyID = localAssemblyContigName.split(":")[0];

        if(assemblyID.length() != 9 )
            throw new IllegalArgumentException("Badly formatted contig name: " + localAssemblyContigName);

        return assemblyID;
    }

    public static SimpleInterval makeOneBpInterval(final String chr, final int pos) {
        Utils.nonNull(chr);
        Utils.validateArg(pos > 0, "given position " + pos + " is non-positive");
        return new SimpleInterval(chr, pos, pos);
    }

    public final Set<String> getSupportingAssemblyIDs() {
        return Collections.unmodifiableSet(
                makeSureAttributeIsList(CONTIG_NAMES)
                        .map(SVLocalContext::parseAssemblyID).collect(Collectors.toSet())
        );
    }
    public final Set<String> getSupporintAssemblyContigNames() {
        return Collections.unmodifiableSet(
                makeSureAttributeIsList(CONTIG_NAMES)
                        .collect(Collectors.toSet())
        );
    }

    /**
     * This exist because for whatever reason,
     * {@link VariantContext#getAttributeAsStringList(String, String)} ()}
     * sometimes returns a giant single string, while
     * {@link VariantContext#getAttributeAsString(String, String)}
     * gives back an array ......
     */
    public final Stream<String> makeSureAttributeIsList(final String attributeKey) {
        return getAttributeAsStringList(attributeKey, "").stream()
                .flatMap(s -> {
                    if ( s.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR) ) {
                        return Arrays.stream( s.split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR) );
                    } else {
                        return Stream.of(s);
                    }
                });
    }

    public abstract boolean isBreakEnd();
    public static boolean isBreakEnd(final VariantContext variantContext) {
        return variantContext
                .getAttributeAsString(SVTYPE, "")
                .equals(BREAKEND_STR);
    }

    abstract Tuple2<SimpleInterval, String> toBedString();


    //==================================================================================================================
    public abstract static class BreakEndSVContext extends SVLocalContext {
        public static final long serialVersionUID = 1L;

        protected BreakEndSVContext(final VariantContext variantContext) {
            super(variantContext);
        }

        @Override
        public final boolean isBreakEnd() {
            return true;
        }

        // TODO: 4/30/18 hack
        public final boolean isUpstreamMate() {
            return getID().endsWith("_1");
        }

        public static Stream<BreakEndSVContext> getOneEnd(final Stream<BreakEndSVContext> breakEndVariants,
                                                                         final boolean firstMate) {
            if (firstMate)
                return breakEndVariants
                        .filter(BreakEndSVContext::isUpstreamMate);// TODO: 3/15/18 insert this key into GATKSVVCFConstant
            else
                return breakEndVariants
                        .filter(bnd -> ! bnd.isUpstreamMate());
        }

        // TODO: 3/16/18 hack for getting the mate location
        public final SimpleInterval getMateRefLoc() {

            final String altSymbAllele = getAlternateAllele(0).toString();

            int i = altSymbAllele.indexOf("]");
            int j;
            if (i >= 0 ) {
                j = altSymbAllele.lastIndexOf("]");
            } else {
                i = altSymbAllele.indexOf("[");
                j = altSymbAllele.lastIndexOf("[");
            }

            return new SimpleInterval(altSymbAllele.substring(i + 1, j));
        }
    }

    /**
     * Exists for {@link #equals(Object)}, {@link #hashCode()} and custom comparator, no other reason.
     * todo: is it possible to talk with {@link SVContext}
     */
    public static final class InvBreakEndContext extends BreakEndSVContext {
        public static final long serialVersionUID = 1L;

        public InvBreakEndContext(final VariantContext variantContext) {
            super(variantContext);
        }

        public boolean isType33() {
            return getAttributeAsBoolean(INV33, false);
        }

        public static Comparator<InvBreakEndContext> makeComparator(final SAMSequenceDictionary refDict) {
            final Comparator<InvBreakEndContext> chrFirst =
                    Comparator.comparingInt(lvc -> refDict.getSequenceIndex(lvc.getContig()));

            return chrFirst
                    .thenComparing(VariantContext::getStart)
                    .thenComparing(VariantContext::getEnd);
        }


        public SimpleInterval getSpanningInterval() {
            int end = getMateRefLoc().getEnd();
            return new SimpleInterval(getContig(), getStart(), end);
        }

        /**
         * Simulating
         * {@link org.broadinstitute.hellbender.utils.test.VariantContextTestUtils#assertVariantContextsAreEqual(VariantContext, VariantContext, List)},
         * but this local version is coded for inversion breakend suspects only,
         * hence testing limited properties and private.
         */
        private static boolean twoVCareEquivalent(final VariantContext v1, final VariantContext v2) {

            if ( !v1.getContig().equals(v2.getContig()) )
                return false;
            if ( v1.getStart() != v2.getStart() )
                return false;
            if ( !v1.getID().equals(v2.getID()) )
                return false;
            if ( !v1.getAlleles().equals(v2.getAlleles()) )
                return false;
            return v1.getAttributes().equals(v2.getAttributes());
        }

        // meant to be overridden
        @Override
        public Tuple2<SimpleInterval, String> toBedString() {
            final String[] splitIdFields = getID().split("_", -1);
            final int start = getStart();
            final Integer end = Integer.valueOf(splitIdFields[splitIdFields.length - 2]);
            final String sz = String.valueOf( end - start);
            final String flag = getAttributeAsBoolean(GATKSVVCFConstants.INV33, false) ? GATKSVVCFConstants.INV33 : GATKSVVCFConstants.INV55;
            final String ctgs = String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR,
                    makeSureAttributeIsList(GATKSVVCFConstants.CONTIG_NAMES).collect(Collectors.toList()));

            final String separator = VCFConstants.INFO_FIELD_SEPARATOR;
            final String annotation = getContig()+("\t")+(start)+("\t")+(end)+("\t")
                    +(flag)+(separator)+(sz)+(separator)+(ctgs);
            return new Tuple2<>( new SimpleInterval(getContig(), start, end), annotation);
        }

        @Override
        public boolean equals(final Object other) {
            if (this == other) return true;
            if (other == null || getClass() != other.getClass()) return false;

            final InvBreakEndContext that = (InvBreakEndContext) other;

            return twoVCareEquivalent(this, that);
        }

        @Override
        public int hashCode() {
            int result = getContig().hashCode();
            result = 31 * result + getStart();
            result = 31 * result + getID().hashCode();
            result = 31 * result + getAlleles().hashCode();
            result = 31 * result + getAttributes().hashCode();
            return result;
        }
    }

    //==================================================================================================================

    public abstract static class SymbolicSVContext extends SVLocalContext {
        public static final long serialVersionUID = 1L;

        SymbolicSVContext(final VariantContext variantContext) {
            super(variantContext);
        }

        @Override
        final public boolean isBreakEnd() {
            return false;
        }
    }
}
