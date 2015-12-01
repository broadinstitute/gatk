package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FastaProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.utils.Utils.calcMD5;

/**
 * Create a dict file from a fasta containing reference sequence.  The output file is a SAM file that contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
@CommandLineProgramProperties(
        summary = "Read fasta or fasta.gz containing reference sequences, and write as a dict file with the sequence dictionary.\n",
        oneLineSummary = "Creates a dict file from reference sequence in fasta format",
        programGroup = FastaProgramGroup.class
)
public final class
CreateSequenceDictionary extends PicardCommandLineProgram {

    // The following attributes define the command-line arguments

    @Argument(
	    shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output a dict file containing only the sequence dictionary")
    public File OUTPUT;

    @Argument(doc = "Put into AS field of sequence dictionary entry if supplied", optional = true)
    public String GENOME_ASSEMBLY;

    @Argument(doc = "Put into UR field of sequence dictionary entry.  If not supplied, input reference file is used",
            optional = true)
    public String URI;

    @Argument(doc = "Put into SP field of sequence dictionary entry", optional = true)
    public String SPECIES;

    @Argument(doc = "Stop after writing this many sequences.  For testing.")
    public int NUM_SEQUENCES = Integer.MAX_VALUE;

    /**
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected String[] customCommandLineValidation() {
        if (REFERENCE_SEQUENCE == null) {
            return new String[]{"********************************\nUSER ERROR: No reference file was specified\n********************************"};
        }
        if (URI == null) {
            URI = "file:" + REFERENCE_SEQUENCE.getAbsolutePath();
        }
        return null;
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     */
    protected Object doWork() {
        if (OUTPUT.exists()) {
            throw new UserException(OUTPUT.getAbsolutePath() + " already exists.  Delete this file and try again, or specify a different output file.");
        }
        final SAMSequenceDictionary sequences = makeSequenceDictionary(REFERENCE_SEQUENCE);
        final SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSequenceDictionary(sequences);
        SAMFileWriter samWriter = null;
        //This writes the header with sequenceDictionary
        try {
            samWriter = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, samHeader, false);
        }
        finally {
            samWriter.close();
        }
        return null;
    }

    /**
     * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     */
    SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, true);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<>();
        final Set<String> sequenceNames = new HashSet<>();
        for (int numSequences = 0; numSequences < NUM_SEQUENCES && (refSeq = refSeqFile.nextSequence()) != null; ++numSequences) {
            if (sequenceNames.contains(refSeq.getName())) {
                throw new UserException.MalformedFile(referenceFile,
                        "Sequence name appears more than once in reference: " + refSeq.getName());
            }
            sequenceNames.add(refSeq.getName());
            ret.add(makeSequenceRecord(refSeq));
        }
        return new SAMSequenceDictionary(ret);
    }

    /**
     * Create one SAMSequenceRecord from a single fasta sequence
     */
    private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
        final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

        // Compute MD5 of upcased bases
        final byte[] bases = refSeq.getBases();
        for (int i = 0; i < bases.length; ++i) {
            bases[i] = StringUtil.toUpperCase(bases[i]);
        }

        ret.setAttribute(SAMSequenceRecord.MD5_TAG, calcMD5(bases));
        if (GENOME_ASSEMBLY != null) {
            ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
        }
        ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
        if (SPECIES != null) {
            ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
        }
        return ret;
    }
}
