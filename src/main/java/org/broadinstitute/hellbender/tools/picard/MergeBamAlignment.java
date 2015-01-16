package org.broadinstitute.hellbender.tools.picard;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.sam.mergealignment.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A command-line tool to merge BAM/SAM alignment info from a third-party aligner with the data in an
 * unmapped BAM file, producing a third BAM file that has alignment data and all the additional data
 * from the unmapped BAM
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Merges alignment data from a SAM or BAM " +
                "file with additional data stored in an unmapped BAM file and produces a third SAM " +
                "or BAM file of aligned and unaligned reads. The purpose is to use information from the " +
                "unmapped BAM to fix up aligner output, so that the resulting file is valid for use by other " +
                "Picard programs. For simple BAM file merges, use MergeSamFiles. NOTE that MergeBamAlignment expects to " +
                "find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it " +
                "to have the same base name as the reference fasta except with the extension '.dict'",
        usageShort = "Merges alignment data from a SAM or BAM with data in an unmapped BAM file",
        programGroup = ReadProgramGroup.class
)
public class MergeBamAlignment extends PicardCommandLineProgram {

    @Deprecated
    @Option(doc = "This argument is ignored and will be removed.", shortName = "PE")
    public Boolean PAIRED_RUN;

    @Option(shortName = "UNMAPPED",
            doc = "Original SAM or BAM file of unmapped reads, which must be in queryname order.")
    public File UNMAPPED_BAM;

    @Option(shortName = "ALIGNED",
            doc = "SAM or BAM file(s) with alignment data.",
            mutex = {"READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM"},
            optional = true)
    public List<File> ALIGNED_BAM;

    @Option(shortName = "R1_ALIGNED",
            doc = "SAM or BAM file(s) with alignment data from the first read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ1_ALIGNED_BAM;

    @Option(shortName = "R2_ALIGNED",
            doc = "SAM or BAM file(s) with alignment data from the second read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ2_ALIGNED_BAM;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Merged SAM or BAM file to write to.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME,
            doc = "Path to the fasta file for the reference sequence.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program group ID of the aligner (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_RECORD_ID;

    @Option(shortName = "PG_VERSION",
            doc = "The version of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

    @Option(shortName = "PG_COMMAND",
            doc = "The command line of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Option(shortName = "PG_NAME",
            doc = "The name of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_NAME;

    @Option(doc = "The expected jump size (required if this is a jumping library). Deprecated. Use EXPECTED_ORIENTATIONS instead",
            shortName = "JUMP",
            mutex = "EXPECTED_ORIENTATIONS",
            optional = true)
    public Integer JUMP_SIZE;

    @Option(doc = "Whether to clip adapters where identified.")
    public boolean CLIP_ADAPTERS = true;

    @Option(doc = "Whether the lane is bisulfite sequence (used when caculating the NM tag).")
    public boolean IS_BISULFITE_SEQUENCE = false;

    @Option(doc = "Whether to output only aligned reads.  ")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc = "The maximum number of insertions or deletions permitted for an alignment to be " +
            "included. Alignments with more than this many insertions or deletions will be ignored. " +
            "Set to -1 to allow any number of insertions or deletions.",
            shortName = "MAX_GAPS")
    public int MAX_INSERTIONS_OR_DELETIONS = 1;

    @Option(doc = "Reserved alignment attributes (tags starting with X, Y, or Z) that should be " +
            "brought over from the alignment data when merging.")
    public List<String> ATTRIBUTES_TO_RETAIN = new ArrayList<String>();

    @Option(doc = "Attributes from the alignment record that should be removed when merging." +
            "  This overrides ATTRIBUTES_TO_RETAIN if they share common tags.")
    public List<String> ATTRIBUTES_TO_REMOVE = new ArrayList<String>();

    @Option(shortName = "R1_TRIM",
            doc = "The number of bases trimmed from the beginning of read 1 prior to alignment")
    public int READ1_TRIM = 0;

    @Option(shortName = "R2_TRIM",
            doc = "The number of bases trimmed from the beginning of read 2 prior to alignment")
    public int READ2_TRIM = 0;

    @Option(shortName = "ORIENTATIONS",
            doc = "The expected orientation of proper read pairs. Replaces JUMP_SIZE",
            mutex = "JUMP_SIZE",
            optional = true)
    public List<SamPairUtil.PairOrientation> EXPECTED_ORIENTATIONS;

    @Option(doc = "Use the aligner's idea of what a proper pair is rather than computing in this program.")
    public boolean ALIGNER_PROPER_PAIR_FLAGS = false;

    @Option(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME,
            doc = "The order in which the merged reads should be output.")
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Option(doc = "Strategy for selecting primary alignment when the aligner has provided more than one alignment " +
            "for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary " +
            "alignment is filtered out for some reason. " +
            "BestMapq expects that multiple alignments will be correlated with HI tag, and prefers the pair of " +
            "alignments with the largest MAPQ, in the absence of a primary selected by the aligner. " +
            "EarliestFragment prefers the alignment which maps the earliest base in the read. Note that EarliestFragment " +
            "may not be used for paired reads. " +
            "BestEndMapq is appropriate for cases in which the aligner is not pair-aware, and does not output the HI tag. " +
            "It simply picks the alignment for each end with the highest MAPQ, and makes those alignments primary, regardless " +
            "of whether the two alignments make sense together." +
            "MostDistant is also for a non-pair-aware aligner, and picks the alignment pair with the largest insert size. " +
            "If all alignments would be chimeric, it picks the alignments for each end with the best MAPQ.  For all algorithms, " +
            "ties are resolved arbitrarily.")
    public PrimaryAlignmentStrategy PRIMARY_ALIGNMENT_STRATEGY = PrimaryAlignmentStrategy.BestMapq;

    @Option(doc = "For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.")
    public boolean CLIP_OVERLAPPING_READS = true;

    @Option(doc = "If false, do not write secondary alignments to output.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = true;

    @Option(shortName = "MC", optional = true, doc = "Adds the mate CIGAR tag (MC) if true, does not if false.")
    public Boolean ADD_MATE_CIGAR = true;

    private static final Log log = Log.getInstance(MergeBamAlignment.class);

    /**
     * Mechanism to bridge between command line option and PrimaryAlignmentSelectionStrategy implementation.
     */
    enum PrimaryAlignmentStrategy {
        BestMapq(BestMapqPrimaryAlignmentSelectionStrategy.class),
        EarliestFragment(EarliestFragmentPrimaryAlignmentSelectionStrategy.class),
        BestEndMapq(BestEndMapqPrimaryAlignmentStrategy.class),
        MostDistant(MostDistantPrimaryAlignmentSelectionStrategy.class);

        private final Class<PrimaryAlignmentSelectionStrategy> clazz;

        PrimaryAlignmentStrategy(final Class<?> clazz) {
            this.clazz = (Class<PrimaryAlignmentSelectionStrategy>) clazz;
        }

        PrimaryAlignmentSelectionStrategy newInstance() {
            try {
                return clazz.newInstance();
            } catch (Exception e) {
                throw new GATKException("Trouble instantiating " + clazz.getName(), e);
            }
        }
    }

    @Override
    protected Object doWork() {
        // Check the files are readable/writable
        SAMProgramRecord prod = null;
        if (PROGRAM_RECORD_ID != null) {
            prod = new SAMProgramRecord(PROGRAM_RECORD_ID);
            prod.setProgramVersion(PROGRAM_GROUP_VERSION);
            prod.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
            prod.setProgramName(PROGRAM_GROUP_NAME);
        }
        // TEMPORARY FIX until internal programs all specify EXPECTED_ORIENTATIONS
        if (JUMP_SIZE != null) {
            EXPECTED_ORIENTATIONS = Arrays.asList(SamPairUtil.PairOrientation.RF);
        } else if (EXPECTED_ORIENTATIONS == null || EXPECTED_ORIENTATIONS.isEmpty()) {
            EXPECTED_ORIENTATIONS = Arrays.asList(SamPairUtil.PairOrientation.FR);
        }

        final SamAlignmentMerger merger = new SamAlignmentMerger(UNMAPPED_BAM, OUTPUT,
                REFERENCE_SEQUENCE, prod, CLIP_ADAPTERS, IS_BISULFITE_SEQUENCE,
                ALIGNED_READS_ONLY, ALIGNED_BAM, MAX_INSERTIONS_OR_DELETIONS,
                ATTRIBUTES_TO_RETAIN, ATTRIBUTES_TO_REMOVE, READ1_TRIM, READ2_TRIM,
                READ1_ALIGNED_BAM, READ2_ALIGNED_BAM, EXPECTED_ORIENTATIONS, SORT_ORDER,
                PRIMARY_ALIGNMENT_STRATEGY.newInstance(), ADD_MATE_CIGAR);
        merger.setClipOverlappingReads(CLIP_OVERLAPPING_READS);
        merger.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
        merger.setKeepAlignerProperPairFlags(ALIGNER_PROPER_PAIR_FLAGS);
        merger.setIncludeSecondaryAlignments(INCLUDE_SECONDARY_ALIGNMENTS);
        merger.mergeAlignment(REFERENCE_SEQUENCE);
        merger.close();

        return null;
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns
     * an array of error messages to be written to the appropriate place.
     */
    protected String[] customCommandLineValidation() {

        if ((PROGRAM_RECORD_ID != null || PROGRAM_GROUP_VERSION != null ||
                PROGRAM_GROUP_COMMAND_LINE != null) &&
                (PROGRAM_RECORD_ID == null || PROGRAM_GROUP_VERSION == null ||
                        PROGRAM_GROUP_COMMAND_LINE == null)) {

            return new String[]{"PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, and " +
                    "PROGRAM_GROUP_COMMAND_LINE must all be supplied or none should " +
                    "be included."};
        }

        final boolean r1sExist = READ1_ALIGNED_BAM != null && READ1_ALIGNED_BAM.size() > 0;
        final boolean r2sExist = READ2_ALIGNED_BAM != null && READ2_ALIGNED_BAM.size() > 0;
        if ((r1sExist && !r2sExist) || (r2sExist && !r1sExist)) {
            return new String[]{"READ1_ALIGNED_BAM and READ2_ALIGNED_BAM " +
                    "must both be supplied or neither should be included.  For " +
                    "single-end read use ALIGNED_BAM."};
        }
        if (ALIGNED_BAM == null || ALIGNED_BAM.size() == 0 && !(r1sExist && r2sExist)) {
            return new String[]{"Either ALIGNED_BAM or the combination of " +
                    "READ1_ALIGNED_BAM and READ2_ALIGNED_BAM must be supplied."};

        }

        return null;
    }

}
