package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils.getCanonicalChromosomes;

public final class TestUtilsForAssemblyBasedSVDiscovery {

    public static final ReferenceMultiSparkSource b37_reference = new ReferenceMultiSparkSource(
            GATKBaseTest.b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b37_seqDict = b37_reference.getReferenceSequenceDictionary(null);
    public static final Set<String> b37_canonicalChromosomes = getCanonicalChromosomes(null, b37_seqDict);
    public static final ReferenceMultiSparkSource b38_reference_chr20_chr21 = new ReferenceMultiSparkSource(
            GATKBaseTest.b38_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b38_seqDict_chr20_chr21 = b38_reference_chr20_chr21.getReferenceSequenceDictionary(null);
    public static final Set<String> b38_canonicalChromosomes = getCanonicalChromosomes(null, b38_seqDict_chr20_chr21);

    /**
     * We are having this because it is SV, especially complex ones, are rare and events on chr20 and 21 are not enough.
     */
    public final static Set<String> hg38CanonicalChromosomes;

    public final static SAMSequenceDictionary bareBoneHg38SAMSeqDict;
    static {
        final List<SAMSequenceRecord> hg38Chromosomes = new ArrayList<>();
        final String hg38ChrBareBoneListFile =  GATKBaseTest.toolsTestDir + "/spark/sv/utils/hg38ChrBareBone.txt";
        try (final Stream<String> records = Files.lines(IOUtils.getPath(( hg38ChrBareBoneListFile )))) {
            records.forEach(line -> {
                final String[] fields = line.split("\t", 2);
                hg38Chromosomes.add(new SAMSequenceRecord(fields[0], Integer.valueOf(fields[1])));
            });
            bareBoneHg38SAMSeqDict = new SAMSequenceDictionary(hg38Chromosomes);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read nonCanonicalContigNamesFile file " + hg38ChrBareBoneListFile, ioe);
        }

        hg38CanonicalChromosomes =
                hg38Chromosomes.subList(0, 24).stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toCollection(LinkedHashSet::new));
    }

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

    public static String makeDummySequence(final char base, final int length) {
        return StringUtils.repeat(base, length);
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

    /**
     * @param primarySAMRecord      input primary line SAM record,
     *                              assumed to have no ambiguity, i.e. no two combinations of its alignments offers equally good story of SV
     */
    public static AssemblyContigWithFineTunedAlignments makeContigAnalysisReady(final String primarySAMRecord,
                                                                                final Set<String> canonicalChromosomes) {
        final AlignedContig alignedContig = fromPrimarySAMRecordString(primarySAMRecord, true);

        final List<AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings> goodAndBadMappings =
                AssemblyContigAlignmentsConfigPicker
                        .pickBestConfigurations(alignedContig, canonicalChromosomes, 0.0);

        return
                AssemblyContigAlignmentsConfigPicker.reConstructContigFromPickedConfiguration(
                        new Tuple2<>(new Tuple2<>(alignedContig.getContigName(), alignedContig.getContigSequence()),
                        goodAndBadMappings))
                .next();
    }
}
