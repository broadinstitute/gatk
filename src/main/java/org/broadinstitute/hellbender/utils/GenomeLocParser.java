package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Factory class for creating GenomeLocs
 */
public final class GenomeLocParser {
    private static final Logger logger = LogManager.getLogger(GenomeLocParser.class);

    public static final String UNMAPPED_LOC_NAME = "unmapped";

    /**
     * How much validation should we do at runtime with this parser?
     */
    public enum ValidationLevel {
        /** Do the standard amount of validation */
        STANDARD,
        /** Don't do any real checking at all */
        NONE
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    // --------------------------------------------------------------------------------------------------------------

    private final MRUCachingSAMSequenceDictionary contigInfo;

    /**
     * How much validation are we doing at runtime with this GenomeLocParser?
     */
    private final ValidationLevel validationLevel;

    /**
     * set our internal reference contig order
     * @param refFile the reference file
     */
    public GenomeLocParser(final ReferenceSequenceFile refFile) {
        this(refFile.getSequenceDictionary());
    }

    /**
     * Create a new GenomeLocParser based on seqDictionary with the standard validation level
     * @param seqDict a non-null sequence dictionary
     */
    public GenomeLocParser(SAMSequenceDictionary seqDict) {
        this(seqDict, ValidationLevel.STANDARD);
    }

    /**
     * Create a genome loc parser based on seqDict with the specified level of validation
     * @param seqDict the sequence dictionary to use when creating genome locs
     * @param validationLevel how much validation should we do of the genome locs at runtime? Purely for testing purposes
     */
    protected GenomeLocParser(SAMSequenceDictionary seqDict, final ValidationLevel validationLevel) {
        Utils.nonNull(validationLevel, "validation level cannot be null");
        if (seqDict == null) { // we couldn't load the reference dictionary
            //logger.info("Failed to load reference dictionary, falling back to lexicographic order for contigs");
            throw new CommandLineException("Failed to load reference dictionary");
        }

        this.validationLevel = validationLevel;
        this.contigInfo = new MRUCachingSAMSequenceDictionary(seqDict);
        if ( logger.isDebugEnabled() ) {
            logger.debug(String.format("Prepared reference sequence contig dictionary"));
            for (SAMSequenceRecord contig : seqDict.getSequences()) {
                logger.debug(String.format(" %s (%d bp)", contig.getSequenceName(), contig.getSequenceLength()));
            }
        }
    }

    /**
     * Determines whether the given contig is valid with respect to the sequence dictionary
     * already installed in the GenomeLoc.
     *
     * @param contig a potentially null string name for the contig
     * @return True if the contig is valid.  False otherwise.
     */
    public final boolean contigIsInDictionary(final String contig) {
        return contig != null && contigInfo.hasContig(contig);
    }

    /**
     * get the contig's SAMSequenceRecord
     *
     * @param contig the string name of the contig
     *
     * @return the sam sequence record
     */
    public final SAMSequenceRecord getContigInfo(final String contig) {
        if ( contig == null || ! contigIsInDictionary(contig) )
            throw new UserException.MalformedGenomeLoc(String.format("Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig));
        return contigInfo.getSequence(contig);
    }

    /**
     * Returns the contig index of a specified string version of the contig
     *
     * @param contig the contig string
     *
     * @return the contig index, -1 if not found
     */
    public final int getContigIndex(final String contig) {
        return getContigInfo(contig).getSequenceIndex();
    }

    protected int getContigIndexWithoutException(final String contig) {
        if ( contig == null || ! contigInfo.hasContig(contig) )
            return -1;
        return contigInfo.getSequenceIndex(contig);
    }

    /**
     * Return the master sequence dictionary used within this GenomeLocParser
     * @return
     */
    public final SAMSequenceDictionary getSequenceDictionary() {
        return contigInfo.getDictionary();
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Low-level creation functions
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * @see #createGenomeLoc(String, int, int, int, boolean) for exact details of the creation.
     *
     * Note that because this function doesn't take the contig index as an argument for contig, it
     * has a slight performance penalty over the version that does take the contig index.  Does not
     * require the created genome loc on the reference genome
     */
    public GenomeLoc createGenomeLoc(String contig, final int start, final int stop) {
        return createGenomeLoc(contig, getContigIndex(contig), start, stop);
    }

    /**
     * @see #createGenomeLoc(String, int, int, int, boolean) for exact details of the creation.
     *
     * Note that because this function doesn't take the contig index as an argument for contig, it
     * has a slight performance penalty over the version that does take the contig index.
     */
    public GenomeLoc createGenomeLoc(final String contig, final int start, final int stop, boolean mustBeOnReference) {
        return createGenomeLoc(contig, getContigIndex(contig), start, stop, mustBeOnReference);
    }

    /**
     * @see #createGenomeLoc(String, int, int, int, boolean) for exact details of the creation.
     *
     * Doesn't require the start and stop to be on the genome
     */
    public GenomeLoc createGenomeLoc(String contig, int index, final int start, final int stop) {
        return createGenomeLoc(contig, index, start, stop, false);
    }

    /**
     * Create a GenomeLoc on contig, starting at start and ending (inclusive) at stop.
     *
     * @param contig the contig name
     * @param index the index into the GATK's SAMSequencingDictionary of contig (passed for efficiency to avoid the lookup)
     * @param start the starting position
     * @param stop  the stop position of this loc, inclusive
     * @param mustBeOnReference if true, this factory will throw a UserException.MalformedGenomeLoc if start or stop isn't on the contig
     *
     * @return a non-null GenomeLoc
     */
    public GenomeLoc createGenomeLoc(final String contig, int index, final int start, final int stop, boolean mustBeOnReference) {
        // optimization: by interning the string we ensure that future comparisons use == not the full string comp
        final String interned = validateGenomeLoc(contig, index, start, stop, mustBeOnReference);
        return new GenomeLoc(interned, index, start, stop);
    }

    /**
     * Create a new GenomeLoc, on contig, including the single position pos.
     *
     * Pos is not required to be on the reference
     *
     * @see #createGenomeLoc(String, int, int, int, boolean) for exact details of the creation.
     *
     * @param contig the contig name
     * @param pos    the start and stop of the created genome loc
     *
     * @return a genome loc representing a single base at the specified postion on the contig
     */
    public GenomeLoc createGenomeLoc(final String contig, final int pos) {
        return createGenomeLoc(contig, getContigIndex(contig), pos, pos);
    }

    /**
     * validate a position or interval on the genome as valid
     *
     * Requires that contig exist in the master sequence dictionary, and that contig index be valid as well.  Requires
     * that start <= stop.
     *
     * if mustBeOnReference is true,
     * performs boundary validation for genome loc INTERVALS:
     * start and stop are on contig and start <= stop
     *
     * @param contig the contig name
     * @param start  the start position
     * @param stop   the stop position
     *
     * @return the interned contig name, an optimization that ensures that contig == the string in the sequence dictionary
     */
    protected String validateGenomeLoc(final String contig, final int contigIndex, final int start, final int stop, final boolean mustBeOnReference) {
        if ( validationLevel == ValidationLevel.NONE )
            return contig;
        else {
            if (stop < start)
                vglHelper(String.format("The stop position %d is less than start %d in contig %s", stop, start, contig));

            final SAMSequenceRecord contigInfo = this.contigInfo.getSequence(contig);
            if ( contigInfo.getSequenceIndex() != contigIndex )
                vglHelper(String.format("The contig index %d is bad, doesn't equal the contig index %d of the contig from a string %s",
                        contigIndex, contigInfo.getSequenceIndex(), contig));

            if ( mustBeOnReference ) {
                if (start < 1)
                    vglHelper(String.format("The start position %d is less than 1", start));

                if (stop < 1)
                    vglHelper(String.format("The stop position %d is less than 1", stop));

                final int contigSize = contigInfo.getSequenceLength();
                if (contigSize == SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH) {
                    logger.warn(String.format("The available sequence dictionary does not contain a sequence length for contig (%s). " +
                            "Skipping validation of the genome loc end coordinate (%d).",
                            contig, stop));
                }
                else if (start > contigSize || stop > contigSize) {
                    vglHelper(String.format("The genome loc coordinates %d-%d exceed the contig size (%d)", start, stop, contigSize));
                }
            }

            return contigInfo.getSequenceName();
        }
    }

    /**
     * Would a genome loc created with the given parameters be valid w.r.t. the master sequence dictionary?
     * @param contig the contig we'd use
     * @param start the start position
     * @param stop the stop
     * @param mustBeOnReference should we require the resulting genome loc to be completely on the reference genome?
     * @return true if this would produce a valid genome loc, false otherwise
     */
    public boolean isValidGenomeLoc(String contig, int start, int stop, boolean mustBeOnReference ) {
        try {
            validateGenomeLoc(contig, getContigIndexWithoutException(contig), start, stop, mustBeOnReference);
            return true;
        } catch ( UserException.MalformedGenomeLoc | GATKException e) {
            return false;
        }
    }

    /**
     * @see #isValidGenomeLoc(String, int, int) with mustBeOnReference == true
     */
    public boolean isValidGenomeLoc(String contig, int start, int stop ) {
        return isValidGenomeLoc(contig, start, stop, true);
    }

    private void vglHelper(final String msg) {
        throw new UserException.MalformedGenomeLoc("Parameters to GenomeLocParser are incorrect:" + msg);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Parsing genome locs
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * parse a genome interval, from a location string
     *
     * Performs interval-style validation:
     *
     * contig is valid; start and stop less than the end; start <= stop, and start/stop are on the contig
     * @param str the string to parse
     *
     * @return a GenomeLoc representing the String
     *
     */
    public GenomeLoc parseGenomeLoc(final String str) {
        try {
            if ( isUnmappedGenomeLocString(str) ) {
                return GenomeLoc.UNMAPPED;
            }

            // Get a list of all possible valid intervals for this query by resolving the query against the sequence
            // dictionary. Throw if there is an ambiguity, otherwise get the unique interval from the list.
            final List<SimpleInterval> allResolvedIntervals = IntervalUtils.getResolvedIntervals(str, contigInfo.getDictionary());
            final Locatable locatable = getUnambiguousInterval(str, allResolvedIntervals);

            final String contig = locatable.getContig();
            final int start = locatable.getStart();
            int stop = locatable.getEnd();

            if (!contigIsInDictionary(contig)) {
                throw new UserException.MalformedGenomeLoc("Contig '" + contig + "' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?");
            }

            if (stop == Integer.MAX_VALUE) {
                // lookup the actual stop position!
                stop = getContigInfo(contig).getSequenceLength();
            }

            return createGenomeLoc(contig, getContigIndex(contig), start, stop, true);
        } catch (UserException.MalformedGenomeLoc e) {
            throw e;
        } catch (IllegalArgumentException | UserException e) {
            throw new UserException.MalformedGenomeLoc("Failed to parse Genome Location string: " + str + ": " + e.getMessage(), e);
        }
    }

    /**
     * Given a (possibly empty) list of valid intervals that have been resolved from a query interval against a
     * sequence dictionary by {@link #parseGenomeLoc}, throw if there are 0, or more than 1, valid intervals in the
     * list, otherwise return the single interval to be used.
     * @param intervalQueryString may not be null
     * @param allResolvedIntervals may not be null
     * @return A single, valid, unambiguous interval
     * @throws UserException.MalformedGenomeLoc if the interval cannot be resolved, or is ambiguous.
     */
    static SimpleInterval getUnambiguousInterval(
            final String intervalQueryString,
            final List<SimpleInterval> allResolvedIntervals)
    {
        Utils.nonNull(intervalQueryString);
        Utils.nonNull(allResolvedIntervals);

        if (allResolvedIntervals.isEmpty()) {
            throw new UserException.MalformedGenomeLoc(
                    String.format("Query interval \"%s\" is not valid for this input.", intervalQueryString));
        } else if (allResolvedIntervals.size() > 1) {
            throw new UserException.MalformedGenomeLoc(
                    String.format(
                            "The query interval \"%s\" is ambiguous and can be interpreted as a query against more than one contig: \"%s\". " +
                                    "The ambiguity can be resolved by providing the interval in a BED file (using zero-based, half-open coordinates).",
                            intervalQueryString,
                            allResolvedIntervals.stream().map(f -> f.getContig()).collect(Collectors.joining(" or "))
                    ));
        } else {
            return allResolvedIntervals.get(0);
        }
    }

    /**
     * @param str String to check
     * @return true if str equals {@link #UNMAPPED_LOC_NAME} (ignoring case)
     */
    public static boolean isUnmappedGenomeLocString(final String str) {
        return str != null && str.trim().equalsIgnoreCase(UNMAPPED_LOC_NAME);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Parsing string representations
    //
    // --------------------------------------------------------------------------------------------------------------

    public GenomeLoc createGenomeLoc( final GATKRead read ) {
        if ( read.isUnmapped() ) {
            return GenomeLoc.UNMAPPED;
        }
        else {
            // Use Math.max to ensure that end >= start (Picard assigns the end to reads that are entirely within an insertion as start-1)
            final int end = Math.max(read.getEnd(), read.getStart());
            return createGenomeLoc(read.getContig(), getSequenceDictionary().getSequenceIndex(read.getContig()), read.getStart(), end, false);
        }
    }

    /**
     * Creates a GenomeLoc from a Tribble feature
     * @param feature
     * @return
     */
    public GenomeLoc createGenomeLoc(final Feature feature) {
        return createGenomeLoc(feature.getContig(), feature.getStart(), feature.getEnd());
    }

    /**
     * Creates a GenomeLoc from a {@link Locatable}.
     *
     * @param locatable the input locatable.
     * @throws IllegalArgumentException if {@code locatable} is {@code null}.
     * @return never {@code null}.
     */
    public GenomeLoc createGenomeLoc(final Locatable locatable) {
        Utils.nonNull(locatable, "the input locatable cannot be null");
        return createGenomeLoc(locatable.getContig(), locatable.getStart(), locatable.getEnd());
    }

    /**
     * Creates a GenomeLoc than spans the entire contig.
     * @param contigName Name of the contig.
     * @return A locus spanning the entire contig.
     */
    public GenomeLoc createOverEntireContig(final String contigName) {
        SAMSequenceRecord contig = contigInfo.getSequence(contigName);
        return createGenomeLoc(contigName,contig.getSequenceIndex(),1,contig.getSequenceLength(), true);
    }

    /**
     * Creates a loc to the left (starting at the loc start + 1) of maxBasePairs size.
     * @param loc The original loc
     * @param maxBasePairs The maximum number of basePairs
     * @return The contiguous loc of up to maxBasePairs length or null if the loc is already at the start of the contig.
     */
    public GenomeLoc createGenomeLocAtStart(final GenomeLoc loc, final int maxBasePairs) {
        if (GenomeLoc.isUnmapped(loc))
            return null;
        final String contigName = loc.getContig();
        final SAMSequenceRecord contig = contigInfo.getSequence(contigName);
        final int contigIndex = contig.getSequenceIndex();

        int start = loc.getStart() - maxBasePairs;
        int stop = loc.getStart() - 1;

        if (start < 1)
            start = 1;
        if (stop < 1)
            return null;

        return createGenomeLoc(contigName, contigIndex, start, stop, true);
    }

    /**
     * Creates a loc padded in both directions by maxBasePairs size (if possible).
     * @param loc      The original loc
     * @param padding  The number of base pairs to pad on either end
     * @return The contiguous loc of length up to the original length + 2*padding (depending on the start/end of the contig).
     */
    public GenomeLoc createPaddedGenomeLoc(final GenomeLoc loc, final int padding) {
        if (GenomeLoc.isUnmapped(loc) || padding == 0)
            return loc;
        else
            return createGenomeLocOnContig(loc.getContig(), loc.getContigIndex(), loc.getStart() - padding, loc.getStop() + padding);
    }

    /**
     * Creates a loc to the right (starting at the loc stop + 1) of maxBasePairs size.
     * @param loc The original loc
     * @param maxBasePairs The maximum number of basePairs
     * @return The contiguous loc of up to maxBasePairs length or null if the loc is already at the end of the contig.
     */
    public GenomeLoc createGenomeLocAtStop(final GenomeLoc loc, final int maxBasePairs) {
        if (GenomeLoc.isUnmapped(loc))
            return null;
        String contigName = loc.getContig();
        SAMSequenceRecord contig = contigInfo.getSequence(contigName);
        int contigIndex = contig.getSequenceIndex();
        int contigLength = contig.getSequenceLength();

        int start = loc.getStop() + 1;
        int stop = loc.getStop() + maxBasePairs;

        if (start > contigLength)
            return null;
        if (stop > contigLength)
            stop = contigLength;

        return createGenomeLoc(contigName, contigIndex, start, stop, true);
    }

    /**
     * @see #createGenomeLocOnContig(String, int, int, int) with the contig index looked up from contig
     */
    public GenomeLoc createGenomeLocOnContig(final String contig, final int start, final int stop) {
        return createGenomeLocOnContig(contig, getContigIndex(contig), start, stop);
    }

    /**
     * Create a new genome loc, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @return a valid genome loc over contig, or null if a meaningful genome loc cannot be created
     */
    public GenomeLoc createGenomeLocOnContig(final String contig, final int contigIndex, final int start, final int stop) {
        final int contigLength = contigInfo.getSequence(contigIndex).getSequenceLength();
        final int boundedStart = Math.max(1, start);
        final int boundedStop = Math.min(contigLength, stop);

        if ( boundedStart > contigLength || boundedStop < 1 )
            // there's no meaningful way to create this genome loc, as the start and stop are off the contig
            return null;
        else
            return createGenomeLoc(contig, contigIndex, boundedStart, boundedStop);
    }
}
