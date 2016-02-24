package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Easy to use creator of artificial BAM files for testing
 *
 * Allows us to make a stream of reads or an index BAM file with read having the following properties
 *
 * - coming from n samples
 * - of fixed read length and aligned to the genome with M operator
 * - having N reads per alignment start
 * - skipping N bases between each alignment start
 * - starting at a given alignment start
 */
public final class ArtificialBAMBuilder {
    private final IndexedFastaSequenceFile reference;
    private final SAMSequenceDictionary dict;

    final int nReadsPerLocus;
    final int nLoci;

    int skipNLoci = 0;
    int alignmentStart = 1;
    int readLength = 10;
    private final ArrayList<String> samples = new ArrayList<>();
    private List<GATKRead> createdReads = null;

    private List<GATKRead> additionalReads = new LinkedList<>();
    private SAMFileHeader header;

    public ArtificialBAMBuilder(int nReadsPerLocus, int nLoci) {
        this(ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000000).getSequenceDictionary(), nReadsPerLocus, nLoci);
    }

    public ArtificialBAMBuilder(final SAMSequenceDictionary dict, int nReadsPerLocus, int nLoci) {
        Utils.nonNull(dict, "dict");
        Utils.validateArg(nReadsPerLocus > 0, "nReadsPerLocus should be positive but was " + nReadsPerLocus);
        Utils.validateArg(nLoci > 0, "nLoci should be positive but was " + nLoci);
        this.nReadsPerLocus = nReadsPerLocus;
        this.nLoci = nLoci;
        this.reference = null;
        this.dict = dict;
        createAndSetHeader(1);
    }

    public IndexedFastaSequenceFile getReference() {
        return reference;
    }

    public ArtificialBAMBuilder createAndSetHeader(final int nSamples) {
        Utils.validateArg(nSamples > 0, "nSamples should be positive but was " + nSamples);

        createdReads = null;
        this.header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(dict);
        samples.clear();

        for ( int i = 0; i < nSamples; i++ ) {
            final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg" + i);
            final String sample = "sample" + i;
            samples.add(sample);
            rg.setSample(sample);
            rg.setPlatform(NGSPlatform.ILLUMINA.getDefaultPlatform());
            header.addReadGroup(rg);
        }

        return this;
    }

    public List<String> getSamples() {
        return samples;
    }

    /**
     * Create a read stream based on the parameters.  The cigar string for each
     * read will be *M, where * is the length of the read.
     *
     * Useful for testing things like LocusIteratorBystate
     *
     * @return a ordered list of reads
     */
    public List<GATKRead> makeReads() {
        if ( createdReads == null ) {
            final String baseName = "read";
            final LinkedList<SAMReadGroupRecord> readGroups = new LinkedList<>();
            for ( final SAMReadGroupRecord rg : header.getReadGroups()) {
                readGroups.add(rg);
            }

            final List<GATKRead> reads = new ArrayList<>(nReadsPerLocus*nLoci);
            for ( int locusI = 0; locusI < nLoci; locusI++) {
                final int locus = locusI * (skipNLoci + 1);
                for ( int readI = 0; readI < nReadsPerLocus; readI++ ) {
                    for ( final SAMReadGroupRecord rg : readGroups ) {
                        final String readName = String.format("%s.%d.%d.%s", baseName, locus, readI, rg.getId());
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, readName, 0, alignmentStart + locus, readLength);
                        read.setReadGroup(rg.getId());
                        reads.add(read);
                    }
                }
            }

            if ( ! additionalReads.isEmpty() ) {
                reads.addAll(additionalReads);
                Collections.sort(reads, new ReadCoordinateComparator(header));
            }

            createdReads = new ArrayList<>(reads);
        }

        return createdReads;
    }

    /**
     * Make an indexed BAM file contains the reads in the builder, marking it for deleteOnExit()
     * @return the BAM file
     */
    public File makeTemporaryBAMFile() throws IOException{
        final File file = IOUtils.createTempFile("tempBAM", ".bam");
        return makeBAMFile(file);
    }

    /**
     * Write the reads from this builder to output, creating an index as well
     * @param output the output BAM file we want to use
     */
    public File makeBAMFile(final File output) {
        Utils.nonNull(output);
        try(final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(output, null, header, false, true, false)){
            for (final GATKRead read : makeReads()) {
                writer.addAlignment(read.convertToSAMRecord(header));
            }
        }
        return output;
    }

    public int getnReadsPerLocus() { return nReadsPerLocus; }

    public int getnLoci() { return nLoci; }

    public int getSkipNLoci() { return skipNLoci; }

    public ArtificialBAMBuilder setSkipNLoci(final int skipNLoci) {
        Utils.validateArg(skipNLoci >= 0, "skipNLoci should be non-negative but was " + skipNLoci);
        this.skipNLoci = skipNLoci;
        createdReads = null;
        return this;
    }

    public int getAlignmentStart() { return alignmentStart; }

    public ArtificialBAMBuilder setAlignmentStart(final int alignmentStart) {
        Utils.validateArg(alignmentStart > 0, "alignmentStart should be positive but was " + alignmentStart);
        this.alignmentStart = alignmentStart;
        createdReads = null;
        return this;
    }

    public int getReadLength() { return readLength; }

    public ArtificialBAMBuilder setReadLength(final int readLength) {
        Utils.validateArg(readLength > 0, "readLength should be positive but was " + readLength);
        this.readLength = readLength;
        createdReads = null;
        return this;
    }

    public SAMFileHeader getHeader() { return header; }

    public ArtificialBAMBuilder setHeader(final SAMFileHeader header) {
        Utils.nonNull(header);
        this.header = header;
        createdReads = null;
        return this;
    }

    public int getAlignmentEnd() {
        return alignmentStart + nLoci * (skipNLoci + 1) + readLength;
    }

    public int getNSamples() { return samples.size(); }

    public int expectedNumberOfReads() {
        return nLoci * nReadsPerLocus * header.getReadGroups().size();
    }

    @Override
    public String toString() {
        return "ArtificialBAMBuilder{" +
                "samples=" + samples +
                ", readLength=" + readLength +
                ", alignmentStart=" + alignmentStart +
                ", skipNLoci=" + skipNLoci +
                ", nLoci=" + nLoci +
                ", nReadsPerLocus=" + nReadsPerLocus +
                '}';
    }
}
