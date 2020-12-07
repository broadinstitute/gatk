package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

import htsjdk.samtools.util.Locatable;
import java.nio.file.Path;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * A BAMWriter that aligns reads to haplotypes and emits their best alignments to a destination.
 *
 */
public class HaplotypeBAMWriter implements AutoCloseable {
    /**
     * Allows us to write out unique names for our synthetic haplotype reads
     */
    private long uniqueNameCounter = 1;

    public static final String DEFAULT_HAPLOTYPE_READ_GROUP_ID = "ArtificialHaplotypeRG";
    // For read filters that need backwards compatibility with the GATK3 artificial haplotype
    public static final String DEFAULT_GATK3_HAPLOTYPE_READ_GROUP_ID = "ArtificialHaplotype";
    private static final int bestHaplotypeMQ = 60;
    private static final int otherMQ = 0;

    private final HaplotypeBAMDestination output;
    private WriterType writerType;
    private boolean writeHaplotypes = true;
    /**
     * Possible modes for writing haplotypes to BAMs
     */
    public enum WriterType {
        /**
         * A mode that's for method developers.  Writes out all of the possible
         * haplotypes considered, as well as reads aligned to each
         */
        ALL_POSSIBLE_HAPLOTYPES,

        /**
         * A mode for users.  Writes out the reads aligned only to the called
         * haplotypes.  Useful to understand why the caller is calling what it is
         */
        CALLED_HAPLOTYPES,

        /**
         * With this option, haplotypes will not be included in the output bam.
         */
        NO_HAPLOTYPES

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
    public HaplotypeBAMWriter(
            final WriterType type,
            final Path outputPath,
            final boolean createBamOutIndex,
            final boolean createBamOutMD5,
            final SAMFileHeader sourceHeader) {

        this(type, new SAMFileDestination(outputPath, createBamOutIndex, createBamOutMD5, sourceHeader, DEFAULT_HAPLOTYPE_READ_GROUP_ID));
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
    public HaplotypeBAMWriter(
            final WriterType type,
            final HaplotypeBAMDestination destination) {

        Utils.nonNull(type, "type cannot be null");
        Utils.nonNull(destination, "destination cannot be null");
        this.writerType = type;
        this.output = destination;
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
     * writerType (ALL_POSSIBLE_HAPLOTYPES or CALLED_HAPLOTYPES) determines inputs to writeHaplotypesAsReads
     * Write out a BAM representing the haplotypes at this site, based on the value for writerType used when
     * the writer was constructed (ALL_POSSIBLE_HAPLOTYPES or CALLED_HAPLOTYPES).
     *
     * @param haplotypes a list of all possible haplotypes at this loc, cannot be null
     * @param paddedReferenceLoc the span of the based reference here, cannot be null
     * @param bestHaplotypes a list of the best (a subset of all) haplotypes that actually went forward into genotyping, cannot be null
     * @param calledHaplotypes a list of the haplotypes that where actually called as non-reference, cannot be null
     * @param readLikelihoods a map from sample -> likelihoods for each read for each of the best haplotypes, cannot be null
     * @param callableRegion the region over which variants are being called
     */
    public void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                              final Locatable paddedReferenceLoc,
                                              final Collection<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                              final Locatable callableRegion) {

        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");
        Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.nonNull(readLikelihoods, "readLikelihoods cannot be null");
        Utils.nonNull(bestHaplotypes, "bestHaplotypes cannot be null");

        if (writerType.equals(WriterType.CALLED_HAPLOTYPES)){
            if (calledHaplotypes.isEmpty()){
                return;
            }
            writeHaplotypesAsReads(calledHaplotypes, calledHaplotypes, paddedReferenceLoc, callableRegion);

        } else if (writerType.equals(WriterType.ALL_POSSIBLE_HAPLOTYPES)){
            writeHaplotypesAsReads(haplotypes, new LinkedHashSet<>(bestHaplotypes), paddedReferenceLoc, callableRegion);
        }

        final int sampleCount = readLikelihoods.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            for (final GATKRead read : readLikelihoods.sampleEvidence(i)) {
                writeReadAgainstHaplotype(read);
            }
        }
    }

    public void writeReadsAlignedToHaplotypes(final Collection<Haplotype> haplotypes,
                                              final Locatable paddedReferenceLoc,
                                              final Collection<Haplotype> bestHaplotypes,
                                              final Set<Haplotype> calledHaplotypes,
                                              final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods) {
        writeReadsAlignedToHaplotypes(haplotypes, paddedReferenceLoc, bestHaplotypes, calledHaplotypes, readLikelihoods, null);
    }

    /**
     * Write out read aligned to haplotype to the BAM file
     *
     * @param read the read we want to write aligned to the reference genome, must not be null
     */
    private void writeReadAgainstHaplotype(final GATKRead read) {
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
     * @param callableRegion the region over which variants are being called
     */
    private void writeHaplotypesAsReads(final Collection<Haplotype> haplotypes,
                                          final Set<Haplotype> bestHaplotypes,
                                          final Locatable paddedReferenceLoc,
                                          final Locatable callableRegion) {
        Utils.nonNull(haplotypes, "haplotypes cannot be null");
        Utils.nonNull(bestHaplotypes, "bestHaplotypes cannot be null");
        Utils.nonNull(paddedReferenceLoc, "paddedReferenceLoc cannot be null");

        if (writeHaplotypes) {
            for (final Haplotype haplotype : haplotypes) {
                writeHaplotype(haplotype, paddedReferenceLoc, bestHaplotypes.contains(haplotype), callableRegion);
            }
        }
    }

    /**
     * Write out a representation of this haplotype as a read
     *
     * @param haplotype a haplotype to write out, must not be null
     * @param paddedRefLoc the reference location, must not be null
     * @param isAmongBestHaplotypes true if among the best haplotypes, false if it was just one possible haplotype
     * @param callableRegion the region over which variants are being called
     */
    private void writeHaplotype(final Haplotype haplotype,
                                final Locatable paddedRefLoc,
                                final boolean isAmongBestHaplotypes,
                                final Locatable callableRegion) {
        Utils.nonNull(haplotype, "haplotype cannot be null");
        Utils.nonNull(paddedRefLoc, "paddedRefLoc cannot be null");

        final SAMRecord record = new SAMRecord(output.getBAMOutputHeader());
        record.setReadBases(haplotype.getBases());
        record.setAlignmentStart(paddedRefLoc.getStart() + haplotype.getAlignmentStartHapwrtRef());
        // Use a base quality value "!" for it's display value (quality value is not meaningful)
        record.setBaseQualities(Utils.dupBytes((byte) '!', haplotype.getBases().length));
        record.setCigar(haplotype.getCigar());
        record.setMappingQuality(isAmongBestHaplotypes ? bestHaplotypeMQ : otherMQ);
        record.setReadName(output.getHaplotypeSampleTag() + uniqueNameCounter++);
        record.setAttribute(output.getHaplotypeSampleTag(), haplotype.hashCode());
        record.setReadUnmappedFlag(false);
        record.setReferenceIndex(output.getBAMOutputHeader().getSequenceIndex(paddedRefLoc.getContig()));
        record.setAttribute(SAMTag.RG.toString(), output.getHaplotypeReadGroupID());
        record.setFlags(SAMFlag.READ_REVERSE_STRAND.intValue());
        if (callableRegion != null) {
            record.setAttribute(AssemblyBasedCallerUtils.CALLABLE_REGION_TAG, callableRegion.toString());
        }

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
