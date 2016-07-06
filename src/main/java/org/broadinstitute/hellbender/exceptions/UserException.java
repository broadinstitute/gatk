package org.broadinstitute.hellbender.exceptions;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.List;
import java.util.Set;

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
        super("A USER ERROR has occurred: " + msg);
    }

    public UserException(final String message, final Throwable throwable) {
        super("A USER ERROR has occurred: " + message, throwable);
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
            super(String.format("Couldn't read file. Error was: %s with exception: %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s. Error was: %s", file.getAbsolutePath(), message));
        }

        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Couldn't read file %s. Error was: %s", file, message));
        }

        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s. Error was: %s with exception: %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e));
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }

    public static class MissingReference extends UserException {
        private static final long serialVersionUID = 0L;

        public MissingReference(String message) { super(message); }
    }

    /**
     * If this exception is thrown during commandline parsing, help will be displayed.
     */
    public static class CommandLineException extends UserException {
        private static final long serialVersionUID = 0L;

        public CommandLineException() {
            super();
        }

        public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
        }
    }

    public static class BadArgumentValue extends CommandLineException {
        private static final long serialVersionUID = 0L;

        public BadArgumentValue(String arg, String value) {
            super(String.format("Argument %s has a bad value: %s", arg, value));
        }

        public BadArgumentValue(String arg, String value, String message){
            super(String.format("Argument %s has a bad value: %s. %s", arg, value,message));
        }

        public BadArgumentValue(String message) {
            super(String.format("Illegal argument value: %s", message));
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
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Couldn't write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), getMessage(e)));
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

        public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }


    public static class BadTmpDir extends UserException {
        private static final long serialVersionUID = 0L;

        public BadTmpDir(String message) {
            super(String.format("Failure working with the tmp directory %s. Override with -Djava.io.tmpdir=X on the command line to a bigger/better file system.  Exact error was %s", System.getProperties().get("java.io.tmpdir"), message));
        }
    }

    public static class MalformedGenomeLoc extends UserException {
        private static final long serialVersionUID = 0L;

        public MalformedGenomeLoc(String message) {
            super(String.format("Badly formed genome loc: %s", message));
        }

        public MalformedGenomeLoc(String message, Exception e) {
            super(String.format("Badly formed genome loc: %s", message), e);
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
            super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)));
        }
    }

    public static class MissingReferenceFaiFile extends UserException {
        private static final long serialVersionUID = 0L;
        public MissingReferenceFaiFile( final File indexFile, final File fastaFile ) {
            super(String.format("Fasta index file %s for reference %s does not exist. Please see %s for help creating it.",
                    indexFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
                    HelpConstants.forumPost("discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference")));
        }
    }

    public static class MissingReferenceDictFile extends UserException {
        private static final long serialVersionUID = 0L;
        public MissingReferenceDictFile( final File dictFile, final File fastaFile ) {
            super(String.format("Fasta dict file %s for reference %s does not exist. Please see %s for help creating it.",
                    dictFile.getAbsolutePath(), fastaFile.getAbsolutePath(),
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

    public static class CannotExecuteRScript extends UserException {
        private static final long serialVersionUID = 0L;
        public CannotExecuteRScript(String message) {
            super(String.format("Unable to execute RScript command: " + message));
        }
    }

    // todo -- fix up exception cause passing
    public static class MissingArgument extends CommandLineException {
        private static final long serialVersionUID = 0L;

        public MissingArgument(String arg, String message) {
            super(String.format("Argument %s was missing: %s", arg, message));
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

        public FailsStrictValidation(File f, ValidateVariants.ValidationType type, String message) {
            super(String.format("File %s fails strict validation: %s of type:", f.getAbsolutePath(), message, type));
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

    public static final class ReferenceAPIReturnedUnexpectedNumberOfBytes extends UserException {
        private static final long serialVersionUID = 0L;
        public ReferenceAPIReturnedUnexpectedNumberOfBytes(final SimpleInterval interval, final byte[] bases) {
            super("Query to genomics service failed for reference interval " + interval + ". Requested " + interval.size() + " bytes but got " + bases.length + ". Perhaps you're querying outside the edge of the contig.");
        }
    }

    public static final class MultipleReferenceSets extends UserException {
        private static final long serialVersionUID = 0L;

        public MultipleReferenceSets(final String referenceSetAssemblyID, final Set<String> referenceSetIds) {
            super("Multiple reference sets found for " + referenceSetAssemblyID + " : " + referenceSetIds + ". Please use a reference set ID that uniquely identifies a reference set.");
        }
    }

    public static final class UnknownReferenceSet extends UserException {
        private static final long serialVersionUID = 0L;

        public UnknownReferenceSet(final String referenceSetAssemblyID) {
            super("There are no known reference set for ID " + referenceSetAssemblyID);
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

        public NoSuitableCodecs(final File file) {
            super("Cannot read " + file + " because no suitable codecs found");
        }
    }

    public static final class WrongFeatureType extends UserException {
        private static final long serialVersionUID = 0L;

        public WrongFeatureType( final File featureFile, final Class<? extends Feature> requiredFeatureType ) {
            super(String.format("File %s is of the wrong type. It should contain Features of type %s, but instead contains Features of a different type.",
                                featureFile.getAbsolutePath(), requiredFeatureType.getSimpleName()));
        }

        public WrongFeatureType( final File featureFile, final Class<? extends Feature> requiredFeatureType, final List<String> actualFeatureTypes ) {
            super(String.format("File %s is of the wrong type. It should contain Features of type %s, but instead contains Features of type(s): %s",
                                featureFile.getAbsolutePath(), requiredFeatureType.getSimpleName(), actualFeatureTypes));
	    }
    }

    public static final class ReadMissingReadGroup extends MalformedBAM {
        private static final long serialVersionUID = 0L;

        public ReadMissingReadGroup(final GATKRead read) {
            super(read, String.format("Read %s is missing the read group (RG) tag, which is required by the GATK.  Please use " + HelpConstants.forumPost("discussion/59/companion-utilities-replacereadgroups to fix this problem"), read.getName()));
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
}
