package org.broadinstitute.hellbender.exceptions;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.List;

/**
 * <p/>
 * Class UserException.
 * <p/>
 * This exception is for errors that are due to user mistakes, such as non-existent or malformed files.
 */
public class UserException extends RuntimeException {
    private static final long serialVersionUID = 0L;

    public UserException() {
        super();
    }

    public UserException(final String msg) {
        super(msg);
    }

    public UserException(final String message, final Throwable throwable) {
        super(message, throwable);
    }

    protected static String getMessage(final Throwable t) {
        final String message = t.getMessage();
        return message != null ? message : t.getClass().getName();
    }

    /**
     * Subtypes of UserException for common kinds of errors
     */

    /**
     * <p/>
     * Class UserException.CouldNotReadInputFile
     * <p/>
     * For generic errors opening/reading from input files
     */
    public static class CouldNotReadInputFile extends UserException {
        private static final long serialVersionUID = 0L;

        public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file. Error was: %s with exception: %s", message, getMessage(e)), e);
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(Path file) {
            super(String.format("Couldn't read file %s", file.toAbsolutePath().toUri()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s. Error was: %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(Path file, String message) {
            super(String.format("Couldn't read file %s. Error was: %s", file.toAbsolutePath().toUri(), message));
        }

        public CouldNotReadInputFile(Path file, String message, Throwable cause) {
            super(String.format("Couldn't read file %s. Error was: %s", file.toAbsolutePath().toUri(), message), cause);
        }

        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Couldn't read file %s. Error was: %s", file, message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s. Error was: %s with exception: %s", file.getAbsolutePath(), message, getMessage(e)), e);
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e), e);
        }

        public CouldNotReadInputFile(Path path, Exception e) {
            this(path, getMessage(e), e);
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }

    public static class MissingReference extends UserException {
        private static final long serialVersionUID = 0L;

        public MissingReference(String message) { super(message); }
    }

    public static class MissingIndex extends UserException {
        private static final long serialVersionUID = 0L;

        public MissingIndex(String file, String message) {
            super(String.format("An index is required but was not found for file %s. %s", file, message));
        }
    }

    public static class CannotHandleGzippedRef extends UserException {
        private static final long serialVersionUID = 0L;

        public CannotHandleGzippedRef() {
            super("The GATK cannot process compressed (.gz) reference sequences. Please unzip the file and try again.  Sorry for the inconvenience.");
        }
    }

    /**
     * <p/>
     * Class UserException.CouldNotCreateOutputFile
     * <p/>
     * For generic errors writing to output files
     */
    public static class CouldNotCreateOutputFile extends UserException {
        private static final long serialVersionUID = 0L;


        public CouldNotCreateOutputFile(File file, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)), e);
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Couldn't write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String file, String message) {
            super(String.format("Couldn't write file %s because %s", file, message));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, getMessage(e)), e);
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), getMessage(e)), e);
        }

        public CouldNotCreateOutputFile(String message, Exception e) {
            super(message, e);
        }
    }

    public static class MalformedRead extends UserException {
        private static final long serialVersionUID = 0L;

        public MalformedRead( GATKRead read, String message ) {
            super(String.format("Read %s is malformed: %s", read, message));
        }
        
        public MalformedRead(String source, String message) {
            super(String.format("Read from source %s is malformed: %s", source, message));
        }
    }

    public static class MisencodedQualityScoresRead extends MalformedRead {
        private static final long serialVersionUID = 0L;

        public MisencodedQualityScoresRead( GATKRead read, String message ) {
            super(read, String.format("Read is using the wrong encoding for quality scores: %s", message));
        }
    }

    public static class BadInput extends UserException {
        private static final long serialVersionUID = 0L;

        public BadInput(String message, Throwable cause){
            super(String.format("Bad input: %s", message), cause);
        }

        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }


    public static class BadTempDir extends UserException {
        private static final long serialVersionUID = 0L;

        private static final String MESSAGE_FORMAT_STRING = "Failure working with the tmp directory %s. Try changing the tmp dir with with --" + StandardArgumentDefinitions.TMP_DIR_NAME + " on the command line.  Exact error was %s";
        public BadTempDir(String message) {
            super(String.format(MESSAGE_FORMAT_STRING, System.getProperties().get("java.io.tmpdir"), message));
        }

        public BadTempDir(String message, Throwable cause) {
            super(String.format(MESSAGE_FORMAT_STRING, System.getProperties().get("java.io.tmpdir"), message), cause);
        }

        public BadTempDir(Path tmpDir, String message, Throwable cause) {
            super(String.format(MESSAGE_FORMAT_STRING, IOUtils.getAbsolutePathWithoutFileProtocol(tmpDir), message), cause);
        }

    }

    public static class MalformedGenomeLoc extends UserException {
        private static final long serialVersionUID = 0L;

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome unclippedLoc: %s", message));
        }

        public MalformedGenomeLoc(String message, Exception e) {
            super(String.format("Badly formed genome unclippedLoc: %s", message), e);
        }
    }

    /**
     * <p/>
     * Class UserException.MalformedFile
     * <p/>
     * For errors parsing files
     */
    public static class MalformedFile extends UserException {
        private static final long serialVersionUID = 0L;

        public MalformedFile(String message) {
            super(String.format("Unknown file is malformed: %s", message));
        }

        public MalformedFile(File f, String message) {
            super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
        }

        public MalformedFile(File f, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)), e);
        }

        public MalformedFile(Path p, String message) {
            super(String.format("File %s is malformed: %s", p.toUri(), message));
        }

        public MalformedFile(Path p, String message, Exception e) {
            super(String.format("File %s is malformed: %s caused by %s", p.toUri(), message, getMessage(e)), e);
        }
    }

    public static class MissingReferenceFaiFile extends UserException {
        private static final long serialVersionUID = 0L;
        public MissingReferenceFaiFile( final Path indexPath, final Path fastaPath ) {
            super(String.format("Fasta index file %s for reference %s does not exist. Please see %s for help creating it.",
                    indexPath.toUri(), fastaPath.toUri(),
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class MissingReferenceDictFile extends UserException {
        private static final long serialVersionUID = 0L;
        public MissingReferenceDictFile( final Path dictFile, final Path fastaFile ) {
            super(String.format("Fasta dict file %s for reference %s does not exist. Please see %s for help creating it.",
                    dictFile.toUri(), fastaFile.toUri(),
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }

        public MissingReferenceDictFile( final String fastaFileName ) {
            super(String.format("Fasta dict file for reference %s does not exist. Please see %s for help creating it.",
                    fastaFileName,
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class IncompatibleRecalibrationTableParameters extends UserException {
        private static final long serialVersionUID = 0L;
        public IncompatibleRecalibrationTableParameters(String s) {
            super(s);
        }
    }

    public static class EmptyIntersection extends UserException {
        private static final long serialVersionUID = 0L;
        public EmptyIntersection(String s) {
            super(s);
        }
    }

    public static class CannotExecuteScript extends UserException {
        private static final long serialVersionUID = 0L;
        public CannotExecuteScript(final String scriptExecutor, String message) {
            super(String.format("Unable to execute %s command: %s", scriptExecutor, message));
        }
    }

    public static class MalformedBAM extends UserException {
        private static final long serialVersionUID = 1l;

        public MalformedBAM(GATKRead read, String message) {
            this("(unknown)", message);
        }

        public MalformedBAM(File file, String message) {
            this(file.toString(), message);
        }

        public MalformedBAM(String source, String message) {
            super(String.format("SAM/BAM/CRAM file %s is malformed: %s", source, message));
        }
    }

    public static class FailsStrictValidation extends UserException {
        private static final long serialVersionUID = 0L;

        public final ValidateVariants.ValidationType type;

        public FailsStrictValidation(String f, ValidateVariants.ValidationType type, String message) {
            super(String.format("Input %s fails strict validation: %s of type:", f, message, type));
            this.type = type;
        }
    }

    public static class MissingContigInSequenceDictionary extends UserException {
        private static final long serialVersionUID = 1L;
        public MissingContigInSequenceDictionary(String contigName, SAMSequenceDictionary dict1) {
            super(String.format("Contig %s not present in the sequence dictionary %s\n",
                    contigName, ReadUtils.prettyPrintSequenceRecords(dict1)));
        }
    }

    /**
     * <p/>
     * Class UserException.MalformedFile
     * <p/>
     * For errors parsing files
     */
    public static class IncompatibleSequenceDictionaries extends UserException {
        private static final long serialVersionUID = 0L;

        public IncompatibleSequenceDictionaries(String message,
                                                String name1,
                                                SAMSequenceDictionary dict1,
                                                String name2,
                                                SAMSequenceDictionary dict2) {
            super(String.format("Input files %s and %s have incompatible contigs: %s.\n  %s contigs = %s\n  %s contigs = %s",
                    name1, name2, message,
                    name1, ReadUtils.prettyPrintSequenceRecords(dict1),
                    name2, ReadUtils.prettyPrintSequenceRecords(dict2)));
        }
    }

    public static class LexicographicallySortedSequenceDictionary extends UserException {
        private static final long serialVersionUID = 0L;

        public LexicographicallySortedSequenceDictionary(String name, SAMSequenceDictionary dict) {
            super(String.format("Lexicographically sorted human genome sequence detected in %s."
                            + "\nFor safety's sake the GATK requires human contigs in karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs."
                            + "\nThis is because all distributed GATK resources are sorted in karyotypic order, and your processing will fail when you need to use these files."
                            + "\nYou can use the ReorderSam utility to fix this problem: " + HelpConstants.forumPost("discussion/58/companion-utilities-reordersam")
                            + "\n  %s contigs = %s",
                    name, name, ReadUtils.prettyPrintSequenceRecords(dict)));
        }
    }

    public static final class Require2BitReferenceForBroadcast extends BadInput {
        private static final long serialVersionUID = 0L;
        public Require2BitReferenceForBroadcast() {
            super("Running this tool with BROADCAST strategy requires a 2bit reference. To create a 2bit reference from an existing fasta file, download faToTwoBit from the link on https://genome.ucsc.edu/goldenPath/help/twoBit.html, then run faToTwoBit in.fasta out.2bit");
	}
    }

    public static final class NoSuitableCodecs extends  UserException {
        private static final long serialVersionUID = 0L;

        public NoSuitableCodecs(final Path file) {
            super("Cannot read " + file + " because no suitable codecs found");
        }
    }

    public static final class WrongFeatureType extends UserException {
        private static final long serialVersionUID = 0L;

        public WrongFeatureType( final File featureFile, final Class<? extends Feature> requiredFeatureType ) {
            super(String.format("File %s is of the wrong type. It should contain Features of type %s, but instead contains Features of a different type.",
                                featureFile.getAbsolutePath(), requiredFeatureType.getSimpleName()));
        }

        public WrongFeatureType( final Path featureFile, final Class<? extends Feature> requiredFeatureType, final List<String> actualFeatureTypes ) {
            super(String.format("File %s is of the wrong type. It should contain Features of type %s, but instead contains Features of type(s): %s",
                    featureFile.toAbsolutePath().toUri(), requiredFeatureType.getSimpleName(), actualFeatureTypes));
        }
    }

    public static final class ReadMissingReadGroup extends MalformedBAM {
        private static final long serialVersionUID = 0L;

        public ReadMissingReadGroup(final GATKRead read) {
            super(read, String.format("Read %s is missing the read group (RG) tag, which is required by the GATK.  Please use " + HelpConstants.forumPost("discussion/59/companion-utilities-replacereadgroups to fix this problem"), read.getName()));
        }
    }

    public static final class HeaderMissingReadGroup extends MalformedBAM {
        private static final long serialVersionUID = 0L;

        public HeaderMissingReadGroup(final GATKRead read) {
            super(read, String.format("Read %s contains an (RG) tag with the group %s which is not found in the file header.", read.getName(), read.getAttributeAsString("RG")));
        }
    }

    public static final class HardwareFeatureException extends UserException {
        private static final long serialVersionUID = 0L;

        public HardwareFeatureException(String message) {
            super(message);
        }

        public HardwareFeatureException(String message, Exception e){
            super(message, e);
        }
    }

    public static final class CouldNotIndexFile extends UserException {
        private static final long serialVersionUID = 0L;

        public CouldNotIndexFile(final File file, final Exception e) {
            super(String.format("Error while trying to create index for %s. Error was: %s: %s",
                    file.getAbsolutePath(), e.getClass().getCanonicalName(), e.getMessage()), e);
        }
    }

    public static final class UnimplementedFeature extends UserException {
        private static final long serialVersionUID = 0L;

        public UnimplementedFeature(String message){
            super(message);
        }
    }
}
