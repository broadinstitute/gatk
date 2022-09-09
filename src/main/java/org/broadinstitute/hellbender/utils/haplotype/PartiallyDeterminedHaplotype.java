package org.broadinstitute.hellbender.utils.haplotype;

import com.google.common.annotations.VisibleForTesting;
import com.netflix.servo.util.Objects;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public final class PartiallyDeterminedHaplotype extends Haplotype {
    private static final long serialVersionUID = 1L;

    //THIS IS DIRECTLY FROM DRAGEN CODE
    final static byte DEL_START = 1;
    final static byte DEL_END = 2;
    final static byte A = 4;
    final static byte C = 8;
    final static byte G = 16;
    final static byte T = 32;

    public byte[] getAlternateBases() {
        return alternateBases;
    }

    private byte[] alternateBases;
    private List<VariantContext> constituentEvents;
    private VariantContext importantVariant;

    /**
     * Compares two haplotypes first by their lengths and then by lexicographic order of their bases.
     */
    public static final Comparator<Haplotype> SIZE_AND_BASE_ORDER =
            Comparator.comparingInt((Haplotype hap) -> hap.getBases().length)
                    .thenComparing(Allele::getBaseString);

    /**
     *
     * @param base
     * @param isRef
     * @param pdBytes
     * @param constituentEvents
     */
    public PartiallyDeterminedHaplotype(final Haplotype base, boolean isRef, byte[] pdBytes, List<VariantContext> constituentEvents, VariantContext eventWithVariant, Cigar cigar) {
        super(base.getBases(), isRef, base.getAlignmentStartHapwrtRef(), cigar);
        this.setGenomeLocation(base.getGenomeLocation());
        this.alternateBases = pdBytes;
        this.constituentEvents = constituentEvents;
        this.importantVariant = eventWithVariant;
    }


    @VisibleForTesting
    public static byte[] getPDBytesForHaplotypes(Allele refAllele, List<Allele> altAlleles) {
        // Asserting we are either indel OR SNP
        byte[] output = new byte[refAllele.length()];
        // SNP case
        if (output.length==1){
            for (Allele a : altAlleles ) {
                switch (a.getBases()[0]) {
                    case 'A':
                        output[0] += A;
                        break;
                    case 'C':
                        output[0] += C;
                        break;
                    case 'T':
                        output[0] += T;
                        break;
                    case 'G':
                        output[0] += G;
                        break;
                    default: throw new RuntimeException("Found unexpected base in alt alleles");
                }
            }
        // DEL case
        } else {
            output[0] += DEL_START;
            output[output.length-1] += DEL_END; //TODO guardrail this
        }
        return output;
    }


    @Override
    public byte[] getDisplayBases() {
        byte[] bytes = getBases().clone();
        for (int i = 0; i < bytes.length; i++) {
            if ((alternateBases[i] & ~3) != 0) { // Don't change the byte if
                bytes[i] = 'P';
            }
        }
        return bytes;
    }

    @Override
    public String toString() {
        String output = "HapLen:"+length() +", "+new String(getDisplayBases());
        output = output + "\nUnresolved Bases["+alternateBases.length+"] "+Arrays.toString(alternateBases);
        return output + "\n"+getCigar().toString()+" "+ constituentEvents.stream()
                .map(v ->(v==this.importantVariant?"*":"")+getDRAGENDebugVariantContextString((int)getStartPosition()).apply(v) )
                .collect(Collectors.joining("->"));
    }

    public static Function<VariantContext, String> getDRAGENDebugVariantContextString(final int offset) {
        return v -> "(" + Integer.toString(v.getStart() - offset) + ",Rlen=" + v.getLengthOnReference() + "," + v.getAlternateAlleles() + ")";
    }

//    /**
//     * Create a new Haplotype derived from this one that exactly spans the provided location
//     *
//     * Note that this haplotype must have a contain a genome loc for this operation to be successful.  If no
//     * GenomeLoc is contained than @throws an IllegalStateException
//     *
//     * Also loc must be fully contained within this Haplotype's genomeLoc.  If not an IllegalArgumentException is
//     * thrown.
//     *
//     * @param loc a location completely contained within this Haplotype's location
//     * @param ignoreRefState should the reference state of the original Haplotype be ignored
//     * @return a new Haplotype within only the bases spanning the provided location, or null for some reason the haplotype would be malformed if
//     */
//    //TODO override this trim
//    public Haplotype trim(final Locatable loc, boolean ignoreRefState) {
//        Utils.nonNull( loc, "Loc cannot be null");
//        Utils.nonNull(genomeLocation, "Cannot trim a Haplotype without containing GenomeLoc");
//        Utils.validateArg(new SimpleInterval(genomeLocation).contains(loc), () -> "Can only trim a Haplotype to a containing span.  My loc is " + genomeLocation + " but wanted trim to " + loc);
//        Utils.nonNull( getCigar(), "Cannot trim haplotype without a cigar " + this);
//
//        final int newStart = loc.getStart() - this.genomeLocation.getStart();
//        final int newStop = newStart + loc.getEnd() - loc.getStart();
//
//        // note: the following returns null if the bases covering the ref interval start or end in a deletion.
//        final byte[] newBases = AlignmentUtils.getBasesCoveringRefInterval(newStart, newStop, getBases(), 0, getCigar());
//
//        if ( newBases == null || newBases.length == 0 ) { // we cannot meaningfully chop down the haplotype, so return null
//            return null;
//        }
//
//        // note: trimCigarByReference does not remove leading or trailing indels, while getBasesCoveringRefInterval does remove bases
//        // of leading and trailing insertions.  We must remove leading and trailing insertions from the Cigar manually.
//        // we keep leading and trailing deletions because these are necessary for haplotypes to maintain consistent reference coordinates
//        final Cigar newCigar = AlignmentUtils.trimCigarByReference(getCigar(), newStart, newStop).getCigar();
//        final boolean leadingInsertion = !newCigar.getFirstCigarElement().getOperator().consumesReferenceBases();
//        final boolean trailingInsertion = !newCigar.getLastCigarElement().getOperator().consumesReferenceBases();
//        final int firstIndexToKeepInclusive = leadingInsertion ? 1 : 0;
//        final int lastIndexToKeepExclusive = newCigar.numCigarElements() - (trailingInsertion ? 1 : 0);
//
//        if (lastIndexToKeepExclusive <= firstIndexToKeepInclusive) {    // edge case of entire cigar is insertion
//            return null;
//        }
//
//        final Cigar leadingIndelTrimmedNewCigar = !(leadingInsertion || trailingInsertion)  ? newCigar :
//                new CigarBuilder(false).addAll(newCigar.getCigarElements().subList(firstIndexToKeepInclusive, lastIndexToKeepExclusive)).make();
//
//        final Haplotype ret = new Haplotype(newBases, ignoreRefState ? false : isReference());
//        ret.setCigar(leadingIndelTrimmedNewCigar);
//        ret.setGenomeLocation(loc);
//        ret.setScore(score);
//        ret.setKmerSize(kmerSize);
//        ret.setAlignmentStartHapwrtRef(newStart + getAlignmentStartHapwrtRef());
//        return ret;
//    }

    @Override
    public boolean equals( final Object h ) {
        return h instanceof PartiallyDeterminedHaplotype
                && getUniquenessValue() == ((Haplotype) h).getUniquenessValue()
                && isReference() == ((Haplotype) h).isReference()
                && Arrays.equals(getBases(), ((Haplotype) h).getBases())
                && Arrays.equals(alternateBases, ((PartiallyDeterminedHaplotype) h).alternateBases);
    }

    @Override
    public int hashCode() {
        // TODO is this good enough for a hash function?
        return Objects.hash(Arrays.hashCode(getBases()),Arrays.hashCode(alternateBases));
    }

}
