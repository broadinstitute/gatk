package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Utility class for test involving reads.
 */
public final class ReadTestUtils {

    /**
     * Returns a list of unpaired reads perfectly matched against a contiguous random location in the
     * reference (non-chimeric reads).
     * <p>
     *     Output reads might map to the forward or reverse strand with equal provability and
     *     they are totally contain in a reference contig (no overhangs).
     * </p>
     * <p>
     *     The resulting reads won't have read-group.
     * </p>
     * <p>
     *     Their names will be result of applying a serial number after the prefix provided.
     * </p>
     * <p>
     *     The resulting read record will be mapped to their original location with the corresponding
     *     single match operation cigar, and their flag values will be those of an unpaired, primary alignment record.
     * </p>
     * <p>
     *     The mapping quality is set to the unknown value {@link SAMRecord#UNKNOWN_MAPPING_QUALITY},
     *     as the reads actually have not been aligned. Also the output reads will lack base qualities.
     * </p>
     * @param rdn random number generator to be used, must not be {@code null}.
     * @param dict reference sequence dictionary, it must match the reference contents in {@code reference}.
     * @param reference reference where the read sequences and locations are based on, must not be {@code null}.
     * @param numberOfReads number of reads to produce, must be 0 or greater.
     * @param namePrefix name prefix for all reads, must not be null.
     * @param serialStart name serial start value, must be 0 or greater.
     * @param minLength minimum read length, must be 1 or greater.
     * @param maxLength maximum read length, must be the same as {@code minLength} or greater.
     * @throws IllegalArgumentException if any of the input parameters do not comply with its requirements.
     * @return never {@code null}, and contains no {@code null} values. It can be further modified by the calling code.
     */
    public static List<SAMRecord> randomErrorFreeUnpairedReads(final Random rdn, final SAMSequenceDictionary dict, final IndexedFastaSequenceFile reference,
                                                             final String namePrefix, final int serialStart,
                                                             final int numberOfReads, final int minLength, final int maxLength) {
        Utils.nonNull(rdn);
        Utils.nonNull(reference);
        Utils.nonNull(namePrefix);
        Utils.nonNull(dict);
        ParamUtils.isPositiveOrZero(serialStart, "serial start");
        ParamUtils.isPositiveOrZero(numberOfReads, "number of reads");
        ParamUtils.isPositive(minLength, "min read length");
        ParamUtils.isPositive(maxLength, "max read length");
        ParamUtils.inRange(maxLength, minLength, Integer.MAX_VALUE, "max read length must be greater or equal to min length");
        final List<SAMRecord> result = new ArrayList<>();
        final SAMFileHeader header = new SAMFileHeader(dict);
        final int numberOfDigitsInSerial = (int) Math.ceil(Math.log10(numberOfReads + serialStart - 1)) + 1;
        final String nameFormat = namePrefix + "%0" + numberOfDigitsInSerial + 'd';
        for (int i = 0; i < numberOfReads; i++) {
            final int contigIndex = rdn.nextInt(dict.size());
            final int seqLength = rdn.nextInt(maxLength - minLength) + minLength;
            final boolean reverse = rdn.nextBoolean();
            final int offset = rdn.nextInt(dict.getSequence(contigIndex).getSequenceLength() - seqLength);
            final int start = offset + 1;
            final byte[] bases = reference.getSubsequenceAt(dict.getSequence(contigIndex).getSequenceName(), start, start + seqLength - 1).getBases();
            if (reverse) {
                SequenceUtil.reverseComplement(bases);
            }
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(String.format(nameFormat, serialStart + i));
            record.setReadNegativeStrandFlag(reverse);
            record.setReadBases(bases);
            record.setReadUnmappedFlag(false);
            record.setReferenceIndex(contigIndex);
            record.setAlignmentStart(start);
            record.setReadPairedFlag(false);
            record.setCigar(new Cigar(Collections.singletonList(new CigarElement(seqLength, CigarOperator.M))));
            record.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY);
            result.add(record);
        }
        return result;
    }

    /**
     * Reads an entire BAM into memory, returning both its SAMFile and all GATKRead records in
     * the bam. Supports both local files and NIO-supported remote filesystems such as GCS.
     *
     * For unit/integration testing purposes only! Do not call this method from actual tools!
     *
     * @param bamPath path or URI to a bam, as a String
     * @return A Pair with the SAMFileHeader as the first element, and a List of all  from the VCF
     *         as the second element
     */
    public static Pair<SAMFileHeader, List<GATKRead>> readEntireBamIntoMemory(final String bamPath) {
        Utils.nonNull(bamPath);

        try ( final ReadsDataSource bamReader = new ReadsDataSource(Paths.get(bamPath)) ) {
            final SAMFileHeader header = bamReader.getHeader();

            final List<GATKRead> reads = new ArrayList<>();
            for ( final GATKRead read : bamReader ) {
                reads.add(read);
            }

            return Pair.of(header, reads);
        }
    }
}
