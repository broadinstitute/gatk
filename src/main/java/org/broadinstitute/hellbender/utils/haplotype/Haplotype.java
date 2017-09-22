package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.Comparator;
import java.util.function.Consumer;

public class Haplotype extends Allele implements Locatable {
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
    public Haplotype trim(final Locatable loc) {
        Utils.nonNull( loc, "Loc cannot be null");
        Utils.nonNull(genomeLocation, "Cannot trim a Haplotype without containing GenomeLoc");
        Utils.validateArg(new SimpleInterval(genomeLocation).contains(loc), () -> "Can only trim a Haplotype to a containing span.  My loc is " + genomeLocation + " but wanted trim to " + loc);
        Utils.nonNull( getCigar(), "Cannot trim haplotype without a cigar " + this);

        final int newStart = loc.getStart() - this.genomeLocation.getStart();
        final int newStop = newStart + loc.getEnd() - loc.getStart();
        final byte[] newBases = AlignmentUtils.getBasesCoveringRefInterval(newStart, newStop, getBases(), 0, getCigar());
        final Cigar newCigar = AlignmentUtils.trimCigarByReference(getCigar(), newStart, newStop);

        if ( newBases == null || AlignmentUtils.startsOrEndsWithInsertionOrDeletion(newCigar) )
            // we cannot meaningfully chop down the haplotype, so return null
        {
            return null;
        }

        final Haplotype ret = new Haplotype(newBases, isReference());
        ret.setCigar(newCigar);
        ret.setGenomeLocation(loc);
        ret.setScore(score);
        ret.setAlignmentStartHapwrtRef(newStart + getAlignmentStartHapwrtRef());
        return ret;
    }

    @Override
    public boolean equals( final Object h ) {
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

    public void setGenomeLocation(final Locatable genomeLocation) {
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
        Utils.validateArg( padSize >= 0, () -> "padSize must be >= 0 but got " + padSize);
        final Cigar extendedHaplotypeCigar = new Cigar(getCigar().getCigarElements());
        if ( padSize > 0 ) {
            extendedHaplotypeCigar.add(new CigarElement(padSize, CigarOperator.M));
        }
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
        Utils.validateArg(!this.cigar.containsOperator(CigarOperator.H), "a haplotype cigar cannot have hard clips");
        Utils.validateArg( this.cigar.getReadLength() == length(), () -> "Read length " + length() + " not equal to the read length of the cigar " + cigar.getReadLength() + " " + this.cigar);
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
    public void setScore(final double score) {
        this.score = score;
    }

    /**
     * Get the span of this haplotype (may be null)
     * @return a potentially null genome loc
     */
    public Locatable getGenomeLocation() {
        return genomeLocation;
    }

    /**
     * Composes a {@link SAMRecord} that contains the sequence and mapping information in the haplotype.
     *
     * @param header the header for the output sam-record file associated with the returned record.
     * @param name the name of the resulting SAMRecord.
     * @param extra additional conversion code, can be {@code null} indicatting that there is no need for additional
     *              conversion.
     *
     * @throws IllegalArgumentException if {@code header} or {@code name} are {@code null}.
     * @return never {@code null}.
     */
    public SAMRecord toSAMRecord(final SAMFileHeader header, final String name, final Consumer<SAMRecord> extra) {
        Utils.nonNull(header, "header cannot be null");
        Utils.nonNull(name);

        final Locatable loc = this.getLocation();
        final Cigar cigar = this.getCigar();
        final Locatable location = this.getLocation();
        final boolean mapped = cigar != null && !cigar.isEmpty() && location != null && location.getContig() != null;

        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadBases(getBases());
        record.setReadUnmappedFlag(!mapped);
        record.setReadPairedFlag(false);
        record.setReadNegativeStrandFlag(false);
        if (mapped) {
            record.setReferenceName(loc.getContig());
            record.setAlignmentStart(location.getStart());
            record.setCigar(cigar);
        }
        if (extra != null) {
            extra.accept(record);
        }
        return record;
    }

    /**
     * Converts the haplotype into a {@link SAMRecord}.
     * <p>
     * This call is equivalent to {@code {@link #toSAMRecord(SAMFileHeader, String, Consumer) toSAMRecord(header, name, null)}}.
     * </p>
     * @see #toSAMRecord(SAMFileHeader, String, Consumer)
     *
     * @param header the header for the output sam-record file associated with the returned record.
     * @param name the name of the resulting SAMRecord.
     *
     * @throws IllegalArgumentException if {@code header} or {@code name} are {@code null}.
     *
     * @return never {@code null}.
     */
    public SAMRecord toSAMRecord(final SAMFileHeader header, final String name) {
        return toSAMRecord(header, name, null);
    }

    /**
     * Returns part of the base's string.
     * @param from 0-based index of the first base to include.
     * @param to 0-based index of the first base not to include.
     * @return never.
     */
    public String getBaseString(final int from, final int to) {
        final byte[] bases = getBases();
        if (to < from) {
            throw new IllegalArgumentException("the to index cannot be smaller than the from index");
        } else if (from < 0) {
            throw new IllegalArgumentException("the from index cannot less than 0");
        } else if (to > bases.length) {
            throw new IllegalArgumentException("the to index cannot be larger than the length of the haplotype");
        } else {
            final StringBuilder builder = new StringBuilder(to - from);
            for (int i = from; i < to; i++) {
                builder.append((char) bases[i]);
            }
            return builder.toString();
        }
    }

    /**
     * Transforms an {@link GATKRead} instance into a {@code Haplotype}.
     * @param record the input record.
     * @param isRef whether the result haplotype is considered reference ({@code true}) or an alternative haplotype ({@code false}).
     * @throws IllegalArgumentException if {@code record} is {@code null}.
     * @return never {@code null}
     */
    public static Haplotype fromGATKRead(final GATKRead record, final boolean isRef) {
        Utils.nonNull(record);
        final Haplotype result = new Haplotype(record.getBases(), isRef);
        if (!record.isUnmapped()) {
            result.setCigar(record.getCigar());
            result.setGenomeLocation(new SimpleInterval(record.getContig(), record.getStart(), record.getEnd()));
        }
        return result;
    }

    @Override
    public String getContig() {
        return genomeLocation == null ? null : genomeLocation.getContig();
    }

    @Override
    public int getStart() {
        return genomeLocation == null ? -1 : genomeLocation.getStart();
    }

    @Override
    public int getEnd() {
        return genomeLocation == null ? -1 : genomeLocation.getEnd();
    }
}
