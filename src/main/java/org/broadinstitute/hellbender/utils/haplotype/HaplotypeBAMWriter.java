package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

import htsjdk.samtools.util.Locatable;
import java.nio.file.Path;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.io.File;
import java.util.Collection;
import java.util.function.Function;
import java.util.Set;

/**
 * A BAMWriter that aligns reads to haplotypes and emits their best alignments to a destination.
 *
 */
public abstract class HaplotypeBAMWriter implements AutoCloseable {
    /**
     * Allows us to write out unique names for our synthetic haplotype reads
     */
    private long uniqueNameCounter = 1;

    public static final String DEFAULT_HAPLOTYPE_READ_GROUP_ID = "ArtificialHaplotypeRG";
    // For read filters that need backwards compatibility with the GATK3 artificial haplotype
    public static final String DEFAULT_GATK3_HAPLOTYPE_READ_GROUP_ID = "ArtificialHaplotype";
    private static final int bestHaplotypeMQ = 60;
    private static final int otherMQ = 0;

    protected final HaplotypeBAMDestination output;
    private boolean writeHaplotypes = true;

    /**
     * Possible modes for writing haplotypes to BAMs
     */
    public enum WriterType {
        /**
         * A mode that's for method developers.  Writes out all of the possible
         * haplotypes considered, as well as reads aligned to each
         */
        ALL_POSSIBLE_HAPLOTYPES(AllHaplotypeBAMWriter::new),

        /**
         * A mode for users.  Writes out the reads aligned only to the called
         * haplotypes.  Useful to understand why the caller is calling what it is
         */
        CALLED_HAPLOTYPES(CalledHaplotypeBAMWriter::new);

        final private Function<HaplotypeBAMDestination, HaplotypeBAMWriter> factory;

        WriterType(Function<HaplotypeBAMDestination, HaplotypeBAMWriter> factory) {
            this.factory = factory;
        }

        /**
         * Create an instance of the  HaplotypeBAMWriter corresponding to this type.
         */
        public HaplotypeBAMWriter create(HaplotypeBAMDestination destination) {
            Utils.nonNull(destination, "destination cannot be null");
            return factory.apply(destination); }
    }

    /**
     * Create a new HaplotypeBAMWriter of type "type" for writing SAMRecords to an output file
     *
     * @param type the type of the writer we want to create, must not be null
     * @param outputPath the destination, must not be null
     * @param createBamOutIndex true to create an index for the bamout
     * @param createBamOutMD5 true to create an MD5 file for the bamout
     * @param sourceHeader the header of the input BAMs used to make calls, must not be null.
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(
            final WriterType type,
            final Path outputPath,
            final boolean createBamOutIndex,
            final boolean createBamOutMD5,
            final SAMFileHeader sourceHeader) {
        Utils.nonNull(type, "type cannot be null");
        Utils.nonNull(outputPath, "outputPath cannot be null");
        Utils.nonNull(sourceHeader, "sourceHeader cannot be null");

        return create(type, new SAMFileDestination(outputPath, createBamOutIndex, createBamOutMD5, sourceHeader, DEFAULT_HAPLOTYPE_READ_GROUP_ID));
    }

    /**
     * Create a new HaplotypeBAMWriter of type "type" for writing SAMRecords to and output destination
     *
     * @param type the type of the writer we want to create
     * @param destination the destination, must not be null. Note that SAM writer associated with the destination must
     *                    have its presorted bit set to false, as reads may come in out of order during writing.
     *
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(final WriterType type, final HaplotypeBAMDestination destination) {
        Utils.nonNull(type, "type cannot be null");
        Utils.nonNull(destination, "destination cannot be null");

        return type.create(destination);
    }

    /**
     * Create a new HaplotypeBAMWriter writing its output to the given destination
     *
     * @param output the output destination, must not be null. Note that SAM writer associated with the destination must
     *               have its presorted bit set to false, as reads may come in out of order during writing.
     */
    protected HaplotypeBAMWriter(final HaplotypeBAMDestination output) {
        Utils.nonNull(output, "output cannot be null");
        this.output = output;
    }

    /**
     * Close any output resource associated with this writer.
     */
    @Override
    public void close() {
        output.close();
    }

    /**
     * Write out a BAM representing for the haplotype caller at this site
     *
     * @param haplotypes a list of all possible haplotypes at this loc
     * @param paddedReferenceLoc the span of the based reference here
     * @param bestHaplotypes a list of the best (a subset of all) haplotypes that actually went forward into genotyping
     * @param calledHaplotypes a list of the haplotypes that where actually called as non-reference
     * @param readLikelihoods a map from sample -> likelihoods for each read for each of the best haplotypes
     */
    public abstract void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                                       final Locatable paddedReferenceLoc,
                                                       final Collection<Haplotype> bestHaplotypes,
                                                       final Set<Haplotype> calledHaplotypes,
                                                       final ReadLikelihoods<Haplotype> readLikelihoods);

    /**
     * Write out read aligned to haplotype to the BAM file
     *
     * @param read the read we want to write aligned to the reference genome, must not be null
     */
    protected void writeReadAgainstHaplotype(final GATKRead read) {
        Utils.nonNull(read, "read cannot be null");
        output.add(read);
    }

    /**
     * If setWriteHaplotypes has been set to true, causes haplotypes to be written to the output destination as reads,
     * marking specifically those that are among the best haplotypes with a higher mapping quality.
     * @param haplotypes a collection of haplotypes to write to the BAM, must not be null
     * @param bestHaplotypes a subset of haplotypes that contains those that are best "either good or called", must not
     *                       be null
     * @param paddedReferenceLoc the genome loc of the padded reference, must not be null
     */
    protected void writeHaplotypesAsReads(final Collection<Haplotype> haplotypes,
                                          final Set<Haplotype> bestHaplotypes,
                                          final Locatable paddedReferenceLoc) {
        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(bestHaplotypes, "bestHaplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");

        if (writeHaplotypes) {
            for (final Haplotype haplotype : haplotypes) {
                writeHaplotype(haplotype, paddedReferenceLoc, bestHaplotypes.contains(haplotype));
            }
        }
    }

    /**
     * Write out a representation of this haplotype as a read
     *
     * @param haplotype a haplotype to write out, must not be null
     * @param paddedRefLoc the reference location, must not be null
     * @param isAmongBestHaplotypes true if among the best haplotypes, false if it was just one possible haplotype
     */
    private void writeHaplotype(final Haplotype haplotype,
                                final Locatable paddedRefLoc,
                                final boolean isAmongBestHaplotypes) {
        Utils.nonNull(haplotype, "haplotype cannot be null");
        Utils.nonNull(paddedRefLoc, "paddedRefLoc cannot be null");

        final SAMRecord record = new SAMRecord(output.getBAMOutputHeader());
        record.setReadBases(haplotype.getBases());
        record.setAlignmentStart(paddedRefLoc.getStart() + haplotype.getAlignmentStartHapwrtRef());
        // Use a base quality value "!" for it's display value (quality value is not meaningful)
        record.setBaseQualities(Utils.dupBytes((byte) '!', haplotype.getBases().length));
        record.setCigar(AlignmentUtils.consolidateCigar(haplotype.getCigar()));
        record.setMappingQuality(isAmongBestHaplotypes ? bestHaplotypeMQ : otherMQ);
        record.setReadName(output.getHaplotypeSampleTag() + uniqueNameCounter++);
        record.setAttribute(output.getHaplotypeSampleTag(), haplotype.hashCode());
        record.setReadUnmappedFlag(false);
        record.setReferenceIndex(output.getBAMOutputHeader().getSequenceIndex(paddedRefLoc.getContig()));
        record.setAttribute(SAMTag.RG.toString(), output.getHaplotypeReadGroupID());
        record.setFlags(SAMFlag.READ_REVERSE_STRAND.intValue());

        output.add(new SAMRecordToGATKReadAdapter(record));
    }

    /**
     * Set the HaplotypeBAMWriter to write out the haplotypes as reads.
     * @param writeHaplotypes true if haplotypes should be written as reads
     */
    public void setWriteHaplotypes(final boolean writeHaplotypes) {
        this.writeHaplotypes = writeHaplotypes;
    }
}
