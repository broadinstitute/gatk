package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils.getCanonicalChromosomes;

public final class SVDiscoveryTestUtilsAndCommonDataProvider {

    // data block ======================================================================================================
    public static final ReferenceMultiSource b37_reference_20_21 = new ReferenceMultiSource(
            GATKBaseTest.b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b37_seqDict_20_21 = b37_reference_20_21.getReferenceSequenceDictionary(null);
    public static final Set<String> b37_canonicalChromosomes = getCanonicalChromosomes(null, b37_seqDict_20_21);
    public static final ReferenceMultiSource b38_reference_chr20_chr21 = new ReferenceMultiSource(
            GATKBaseTest.b38_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b38_seqDict_chr20_chr21 = b38_reference_chr20_chr21.getReferenceSequenceDictionary(null);
    public static final Set<String> b38_canonicalChromosomes = getCanonicalChromosomes(null, b38_seqDict_chr20_chr21);

    // utils block =====================================================================================================

    public static byte[] getReverseComplimentCopy(final byte[] sequence) {
        final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);
        SequenceUtil.reverseComplement(sequenceCopy);
        return sequenceCopy;
    }

    public static byte[] makeDummySequence(final int length, byte base) {
        final byte[] result = new byte[length];
        Arrays.fill(result, base);
        return result;
    }

    // WARNING: THIS SHOULD BE USED ONLY FOR CONSTRUCTING ALIGNMENT INTERVAL FOR SV TESTS FROM WELL FORMATTED SAM STRING
    public static AlignmentInterval fromSAMRecordString(final String samRecordStringWithExplicitTabEscapeSequenceAndNoXAField, final boolean hasSATag) {
        final String[] fields = samRecordStringWithExplicitTabEscapeSequenceAndNoXAField.split("\t");
        final int samFlag = Integer.valueOf( fields[1] ) ;
        final String chr = fields[2];
        final int start = Integer.valueOf(fields[3]);
        final int mapQual = Integer.valueOf( fields[4] );
        final Cigar cigar = TextCigarCodec.decode(fields[5]);

        final int idx = hasSATag ? 11 : 10;
        final int numMismatch = Integer.valueOf( fields[idx+3].substring(3 + fields[idx+3].indexOf(":i:")) );
        final int alignerScore = Integer.valueOf( fields[idx+4].substring(3 + fields[idx+4].indexOf(":i:")) );

        final boolean forwardStrand = SAMFlag.READ_REVERSE_STRAND.isUnset(samFlag);
        final Cigar readCigar = forwardStrand ? cigar : CigarUtils.invertCigar(cigar);

        final SimpleInterval refSpan = new SimpleInterval(chr, start, start - 1 + cigar.getReferenceLength());
        return new AlignmentInterval(refSpan,
                SvCigarUtils.getNumClippedBases(true, readCigar) + 1,
                SvCigarUtils.getUnclippedReadLength(readCigar) - SvCigarUtils.getNumClippedBases(false, readCigar),
                readCigar,
                forwardStrand,
                mapQual, numMismatch, alignerScore,
                ContigAlignmentsModifier.AlnModType.NONE);
    }

    // WARNING: THIS SHOULD BE USED ONLY FOR CONSTRUCTING ALIGNMENT INTERVAL FOR SV TESTS FROM WELL FORMATTED SAM STRING OF A PRIMARY RECORD THAT HAS NO XA TAG
    public static AlignedContig fromPrimarySAMRecordString(final String samRecordStringWithExplicitTabEscapeSequenceAndNoXAField,
                                                           final boolean hasSATag) {
        final AlignmentInterval primaryAlignment = fromSAMRecordString(samRecordStringWithExplicitTabEscapeSequenceAndNoXAField, hasSATag);

        final String[] fields = samRecordStringWithExplicitTabEscapeSequenceAndNoXAField.split("\t");
        final String readName = fields[0];
        final int samFlag = Integer.valueOf( fields[1] ) ;
        final byte[] sequence = fields[9].getBytes();
        final boolean forwardStrand = SAMFlag.READ_REVERSE_STRAND.isUnset(samFlag);
        if (!forwardStrand) {
            SequenceUtil.reverseComplement(sequence);
        }

        if (hasSATag) {
            final String saTag = fields[11];
            final String[] supplements = saTag.substring(3 + saTag.indexOf(":Z:")).split(";");
            final List<AlignmentInterval> alignments = new ArrayList<>(supplements.length + 1);
            alignments.add(primaryAlignment);
            for (final String sup : supplements) {
                alignments.add( new AlignmentInterval(sup) );
            }
            return new AlignedContig(readName, sequence, alignments);
        } else {
            return new AlignedContig(readName, sequence, Collections.singletonList(primaryAlignment));
        }
    }

}
