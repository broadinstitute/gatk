package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.HasGenomeLocation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

public final class Haplotype extends Allele {
    private static final long serialVersionUID = 1L;

    /**
     * Compares two haplotypes first by their lengths and then by lexicographic order of their bases.
     */
    public static final Comparator<Haplotype> SIZE_AND_BASE_ORDER =
            Comparator.comparingInt((Haplotype hap) -> hap.getBases().length)
                      .thenComparing(hap -> hap.getBaseString());

    private Locatable genomeLocation = null;
    private EventMap eventMap = null;
    private Cigar cigar;
    private int alignmentStartHapwrtRef;
    private double score = Double.NaN;

    /**
     * Main constructor
     *
     * @param bases a non-null array of bases
     * @param isRef is this the reference haplotype?
     */
    public Haplotype( final byte[] bases, final boolean isRef ) {
        super(Arrays.copyOf(bases, bases.length), isRef);
    }

    /**
     * Create a new non-ref haplotype
     *
     * @param bases a non-null array of bases
     */
    public Haplotype( final byte[] bases ) {
        this(bases, false);
    }

    /**
     * Create a new haplotype with bases
     *
     * Requires bases.length == cigar.getReadLength()
     *
     * @param bases a non-null array of bases
     * @param isRef is this the reference haplotype?
     * @param alignmentStartHapwrtRef offset of this haplotype w.r.t. the reference
     * @param cigar the cigar that maps this haplotype to the reference sequence
     */
    public Haplotype( final byte[] bases, final boolean isRef, final int alignmentStartHapwrtRef, final Cigar cigar) {
        this(bases, isRef);
        this.alignmentStartHapwrtRef = alignmentStartHapwrtRef;
        setCigar(cigar);
    }

    public Haplotype( final byte[] bases, final Locatable loc ) {
        this(bases, false);
        this.genomeLocation = loc;
    }

    /**
     * Create a new Haplotype derived from this one that exactly spans the provided location
     *
     * Note that this haplotype must have a contain a genome loc for this operation to be successful.  If no
     * GenomeLoc is contained than @throws an IllegalStateException
     *
     * Also loc must be fully contained within this Haplotype's genomeLoc.  If not an IllegalArgumentException is
     * thrown.
     *
     * @param loc a location completely contained within this Haplotype's location
     * @return a new Haplotype within only the bases spanning the provided location, or null for some reason the haplotype would be malformed if
     */
    public Haplotype trim(final GenomeLoc loc) {
        if ( loc == null ) throw new IllegalArgumentException("Loc cannot be null");
        if ( genomeLocation == null ) throw new IllegalStateException("Cannot trim a Haplotype without containing GenomeLoc");
        if ( ! new SimpleInterval(genomeLocation).contains(loc) ) throw new IllegalArgumentException("Can only trim a Haplotype to a containing span.  My loc is " + genomeLocation + " but wanted trim to " + loc);
        if ( getCigar() == null ) throw new IllegalArgumentException("Cannot trim haplotype without a cigar " + this);

        final int newStart = loc.getStart() - this.genomeLocation.getStart();
        final int newStop = newStart + loc.size() - 1;
        final byte[] newBases = AlignmentUtils.getBasesCoveringRefInterval(newStart, newStop, getBases(), 0, getCigar());
        final Cigar newCigar = AlignmentUtils.trimCigarByReference(getCigar(), newStart, newStop);

        if ( newBases == null || AlignmentUtils.startsOrEndsWithInsertionOrDeletion(newCigar) )
            // we cannot meaningfully chop down the haplotype, so return null
            return null;

        final Haplotype ret = new Haplotype(newBases, isReference());
        ret.setCigar(newCigar);
        ret.setGenomeLocation(loc);
        ret.setAlignmentStartHapwrtRef(newStart + getAlignmentStartHapwrtRef());
        return ret;
    }

    @Override
    public boolean equals( Object h ) {
        return h instanceof Haplotype && Arrays.equals(getBases(), ((Haplotype) h).getBases());
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(getBases());
    }

    public EventMap getEventMap() {
        return eventMap;
    }

    public void setEventMap( final EventMap eventMap ) {
        this.eventMap = eventMap;
    }

    @Override
    public String toString() {
        return getDisplayString();
    }

    /**
     * Get the span of this haplotype (may be null)
     * @return a potentially null genome loc
     */
    public Locatable getLocation() {
        return this.genomeLocation;
    }

    public void setGenomeLocation(GenomeLoc genomeLocation) {
        this.genomeLocation = genomeLocation;
    }

    public long getStartPosition() {
        return genomeLocation.getStart();
    }

    public long getStopPosition() {
        return genomeLocation.getEnd();
    }

    public int getAlignmentStartHapwrtRef() {
        return alignmentStartHapwrtRef;
    }

    public void setAlignmentStartHapwrtRef( final int alignmentStartHapwrtRef ) {
        this.alignmentStartHapwrtRef = alignmentStartHapwrtRef;
    }

    /**
     * Get the cigar for this haplotype.  Note that the cigar is guaranteed to be consolidated
     * in that multiple adjacent equal operates will have been merged
     * @return the cigar of this haplotype
     */
    public Cigar getCigar() {
        return cigar;
    }

    /**
     * Get the haplotype cigar extended by padSize M at the tail, consolidated into a clean cigar
     *
     * @param padSize how many additional Ms should be appended to the end of this cigar.  Must be >= 0
     * @return a newly allocated Cigar that consolidate(getCigar + padSize + M)
     */
    public Cigar getConsolidatedPaddedCigar(final int padSize) {
        if ( padSize < 0 ) throw new IllegalArgumentException("padSize must be >= 0 but got " + padSize);
        final Cigar extendedHaplotypeCigar = new Cigar(getCigar().getCigarElements());
        if ( padSize > 0 ) extendedHaplotypeCigar.add(new CigarElement(padSize, CigarOperator.M));
        return AlignmentUtils.consolidateCigar(extendedHaplotypeCigar);
    }

    /**
     * Set the cigar of this haplotype to cigar.
     *
     * Note that this function consolidates the cigar, so that 1M1M1I1M1M => 2M1I2M
     *
     * @param cigar a cigar whose readLength == length()
     */
    public void setCigar( final Cigar cigar ) {
        this.cigar = AlignmentUtils.consolidateCigar(cigar);
        if ( this.cigar.getReadLength() != length() )
            throw new IllegalArgumentException("Read length " + length() + " not equal to the read length of the cigar " + cigar.getReadLength() + " " + this.cigar);
    }

    public Haplotype insertAllele( final Allele refAllele, final Allele altAllele, final int refInsertLocation, final int genomicInsertLocation ) {
        // refInsertLocation is in ref haplotype offset coordinates NOT genomic coordinates
        final int haplotypeInsertLocation = ReadUtils.getReadCoordinateForReferenceCoordinate(alignmentStartHapwrtRef, cigar, refInsertLocation, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        final byte[] myBases = this.getBases();
        if( haplotypeInsertLocation == -1 || haplotypeInsertLocation + refAllele.length() >= myBases.length ) { // desired change falls inside deletion so don't bother creating a new haplotype
            return null;
        }

        byte[] newHaplotypeBases = {};
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, ArrayUtils.subarray(myBases, 0, haplotypeInsertLocation)); // bases before the variant
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, altAllele.getBases()); // the alt allele of the variant
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, ArrayUtils.subarray(myBases, haplotypeInsertLocation + refAllele.length(), myBases.length)); // bases after the variant
        return new Haplotype(newHaplotypeBases);
    }

    /**
     * Get the score (an estimate of the support) of this haplotype
     * @return a double, where higher values are better
     */
    public double getScore() {
        return score;
    }

    /**
     * Set the score (an estimate of the support) of this haplotype.
     *
     * Note that if this is the reference haplotype it is always given Double.MAX_VALUE score
     *
     * @param score a double, where higher values are better
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**
     * Get the span of this haplotype (may be null)
     * @return a potentially null genome loc
     */
    public Locatable getGenomeLocation() {
        return genomeLocation;
    }

    public static Map<Allele,Haplotype> makeHaplotypeListFromAlleles(final List<Allele> alleleList,
                                                                               final int startPos,
                                                                               final ReferenceContext ref,
                                                                               final int haplotypeSize,
                                                                               final int numPrefBases) {

        LinkedHashMap<Allele,Haplotype> haplotypeMap = new LinkedHashMap<>();

        Allele refAllele = null;

        for (Allele a:alleleList) {
            if (a.isReference()) {
                refAllele = a;
                break;
            }
        }

        if (refAllele == null)
            throw new GATKException("BUG: no ref alleles in input to makeHaplotypeListfrom Alleles at loc: "+ startPos);

        final byte[] refBases = ref.getBases();

        final int startIdxInReference = 1 + startPos - numPrefBases - ref.getWindow().getStart();
        final String basesBeforeVariant = new String(Arrays.copyOfRange(refBases, startIdxInReference, startIdxInReference + numPrefBases));

        // protect against long events that overrun available reference context
        final int startAfter = Math.min(startIdxInReference + numPrefBases + refAllele.getBases().length - 1, refBases.length);
        final String basesAfterVariant = new String(Arrays.copyOfRange(refBases, startAfter, refBases.length));

        // Create location for all haplotypes
        final int startLoc = ref.getWindow().getStart() + startIdxInReference;
        final int stopLoc = startLoc + haplotypeSize-1;

        final Locatable locus = new SimpleInterval(ref.getInterval().getContig(),startLoc,stopLoc);

        for (final Allele a : alleleList) {

            final byte[] alleleBases = a.getBases();
            // use string concatenation
            String haplotypeString = basesBeforeVariant + new String(Arrays.copyOfRange(alleleBases, 1, alleleBases.length)) + basesAfterVariant;
            haplotypeString = haplotypeString.substring(0,haplotypeSize);

            haplotypeMap.put(a,new Haplotype(haplotypeString.getBytes(), locus));
        }

        return haplotypeMap;
    }
}
