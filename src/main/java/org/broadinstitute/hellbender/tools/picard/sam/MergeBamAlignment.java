package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SamPairUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.mergealignment.*;

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
        summary = "Merges alignment data from a SAM/BAM/CRAM " +
                "file with additional data stored in an unmapped SAM/BAM/CRAM file and produces a third SAM/BAM/CRAM " +
                "file of aligned and unaligned reads. The purpose is to use information from the " +
                "unmapped SAM/BAM/CRAM to fix up aligner output, so that the resulting file is valid for use by other " +
                "Picard programs. For simple SAM/BAM/CRAM file merges, use MergeSamFiles. NOTE that MergeBamAlignment expects to " +
                "find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it " +
                "to have the same base name as the reference fasta except with the extension '.dict'",
        oneLineSummary = "Merges alignment data from a SAM/BAM with data in an unmapped SAM/BAM/CRAM file",
        programGroup = ReadProgramGroup.class
)
public final class MergeBamAlignment extends PicardCommandLineProgram {

    @Argument(shortName = "UNMAPPED",
            doc = "Original SAM/BAM/CRAM file of unmapped reads, which must be in queryname order.")
    public File UNMAPPED_BAM;

    @Argument(shortName = "ALIGNED", fullName = "ALIGNED_BAM",
            doc = "SAM/BAM/CRAM file(s) with alignment data.",
            mutex = {"READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM"},
            optional = true)
    public List<File> ALIGNED_BAM;

    @Argument(shortName = "R1_ALIGNED", fullName= "READ1_ALIGNED_BAM",
            doc = "SAM/BAM/CRAM file(s) with alignment data from the first read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ1_ALIGNED_BAM;

    @Argument(shortName = "R2_ALIGNED", fullName= "READ2_ALIGNED_BAM",
            doc = "SAM/BAM file(s) with alignment data from the second read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ2_ALIGNED_BAM;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Merged SAM/BAM/CRAM file to write to.")
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program group ID of the aligner (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_RECORD_ID;

    @Argument(shortName = "PG_VERSION",
            doc = "The version of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

    @Argument(shortName = "PG_COMMAND",
            doc = "The command line of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Argument(shortName = "PG_NAME",
            doc = "The name of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_NAME;

    @Argument(doc = "Whether to clip adapters where identified.")
    public boolean CLIP_ADAPTERS = true;

    @Argument(doc = "Whether the lane is bisulfite sequence (used when caculating the NM tag).")
    public boolean IS_BISULFITE_SEQUENCE = false;

    @Argument(doc = "Whether to output only aligned reads.  ")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc = "The maximum number of insertions or deletions permitted for an alignment to be " +
            "included. Alignments with more than this many insertions or deletions will be ignored. " +
            "Set to -1 to allow any number of insertions or deletions.",
            shortName = "MAX_GAPS")
    public int MAX_INSERTIONS_OR_DELETIONS = 1;

    @Argument(doc = "Reserved alignment attributes (tags starting with X, Y, or Z) that should be " +
            "brought over from the alignment data when merging.",
             optional = true)
    public List<String> ATTRIBUTES_TO_RETAIN = new ArrayList<>();

    @Argument(doc = "Attributes from the alignment record that should be removed when merging." +
            "  This overrides ATTRIBUTES_TO_RETAIN if they share common tags.",
            optional = true)
    public List<String> ATTRIBUTES_TO_REMOVE = new ArrayList<>();

    @Argument(shortName = "R1_TRIM",
            doc = "The number of bases trimmed from the beginning of read 1 prior to alignment")
    public int READ1_TRIM = 0;

    @Argument(shortName = "R2_TRIM",
            doc = "The number of bases trimmed from the beginning of read 2 prior to alignment")
    public int READ2_TRIM = 0;

    @Argument(shortName = "ORIENTATIONS",
            doc = "The expected orientation of proper read pairs.",
            optional = true)
    public List<SamPairUtil.PairOrientation> EXPECTED_ORIENTATIONS;

    @Argument(doc = "Use the aligner's idea of what a proper pair is rather than computing in this program.")
    public boolean ALIGNER_PROPER_PAIR_FLAGS = false;

    @Argument(shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME,
            doc = "The order in which the merged reads should be output.")
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Argument(doc = "Strategy for selecting primary alignment when the aligner has provided more than one alignment " +
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

    @Argument(doc = "For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.")
    public boolean CLIP_OVERLAPPING_READS = true;

    @Argument(doc = "If false, do not write secondary alignments to output.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = true;

    @Argument(shortName = "MC", optional = true, doc = "Adds the mate CIGAR tag (MC) if true, does not if false.")
    public Boolean ADD_MATE_CIGAR = true;

    /**
     * Mechanism to bridge between command line option and PrimaryAlignmentSelectionStrategy implementation.
     */
    enum PrimaryAlignmentStrategy {
        BestMapq(BestMapqPrimaryAlignmentSelectionStrategy.class),
        EarliestFragment(EarliestFragmentPrimaryAlignmentSelectionStrategy.class),
        BestEndMapq(BestEndMapqPrimaryAlignmentStrategy.class),
        MostDistant(MostDistantPrimaryAlignmentSelectionStrategy.class);

        private final Class<? extends PrimaryAlignmentSelectionStrategy> clazz;

        PrimaryAlignmentStrategy(final Class<? extends PrimaryAlignmentSelectionStrategy> clazz) {
            this.clazz = clazz;
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
        if (EXPECTED_ORIENTATIONS == null || EXPECTED_ORIENTATIONS.isEmpty()) {
            EXPECTED_ORIENTATIONS = Arrays.asList(SamPairUtil.PairOrientation.FR);
        }

        final SamAlignmentMerger merger = new SamAlignmentMerger(UNMAPPED_BAM, OUTPUT,
                REFERENCE_SEQUENCE, prod, CLIP_ADAPTERS, IS_BISULFITE_SEQUENCE,
                ALIGNED_READS_ONLY, ALIGNED_BAM, MAX_INSERTIONS_OR_DELETIONS,
                ATTRIBUTES_TO_RETAIN, ATTRIBUTES_TO_REMOVE, READ1_TRIM, READ2_TRIM,
                READ1_ALIGNED_BAM, READ2_ALIGNED_BAM, EXPECTED_ORIENTATIONS, SORT_ORDER,
                PRIMARY_ALIGNMENT_STRATEGY.newInstance(), ADD_MATE_CIGAR);
        merger.setClipOverlappingReads(CLIP_OVERLAPPING_READS);
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
