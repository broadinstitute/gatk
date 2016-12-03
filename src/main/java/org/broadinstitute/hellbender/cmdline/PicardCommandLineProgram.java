package org.broadinstitute.hellbender.cmdline;

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;

/**
 * Base class for all Picard tools. Includes standard options for dealing with common sequence data formats.
 */
public abstract class PicardCommandLineProgram extends CommandLineProgram {

    @Argument(doc = "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT " +
            "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
            "do not otherwise need to be decoded.", common=true)
    public ValidationStringency VALIDATION_STRINGENCY = ValidationStringency.DEFAULT_STRINGENCY;

    @Argument(doc = "Compression level for all compressed files created (e.g. BAM and GELI).", common=true)
    public int COMPRESSION_LEVEL = BlockCompressedStreamConstants.DEFAULT_COMPRESSION_LEVEL;

    @Argument(doc = "When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.", optional=true, common=true)
    public Integer MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();

    @Argument(doc = "Whether to create a BAM index when writing a coordinate-sorted BAM file.", common=true)
    public Boolean CREATE_INDEX = Defaults.CREATE_INDEX;

    @Argument(doc="Whether to create an MD5 digest for any BAM or FASTQ files created.  ", common=true)
    public boolean CREATE_MD5_FILE = Defaults.CREATE_MD5;

    @Argument(fullName =  StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file.", common = true, optional = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    @Override
    public Object instanceMain(final String[] argv) {
        // First, we parse the commandline arguments, then we set important statics like VALIDATION_STRINGENCY, and
        // finally, we call into the normal instance main (post arg-parsing). If don't start with the argument parsing
        // we always get default values for VALIDATION_STRINGENCY, COMPRESSION_LEVEL, etc.
        if (!parseArgs(argv)) {
            //an information only argument like help or version was specified, just exit
            return 0;
        }
        // set general SAM/BAM parameters
        SamReaderFactory.setDefaultValidationStringency(VALIDATION_STRINGENCY);

        BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL);

        if (MAX_RECORDS_IN_RAM != null) {
            SAMFileWriterImpl.setDefaultMaxRecordsInRam(MAX_RECORDS_IN_RAM);
        }

        if (CREATE_INDEX){
            SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
        }

        SAMFileWriterFactory.setDefaultCreateMd5File(CREATE_MD5_FILE);

        // defer to parent to finish the initialization and starting the program.
        return instanceMainPostParseArgs();
    }

    /**
     * Create a common SAMFileWriter for use with Picard tools.
     *
     * @param outputFile    - if this file has a .cram extension then a reference is required. Can not be null.
     * @param referenceFile - the reference source to use. Can not be null if a output file has a .cram extension.
     * @param header        - header to be used for the output writer
     * @param preSorted     - if true then the records must already be sorted to match the header sort order
     * @return SAMFileWriter
     */
    public SAMFileWriter createSAMWriter(
            final File outputFile,
            final File referenceFile,
            final SAMFileHeader header,
            final boolean preSorted)
    {
        BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL);

        SAMFileWriterFactory factory = new SAMFileWriterFactory()
                .setCreateIndex(CREATE_INDEX)
                .setCreateMd5File(CREATE_MD5_FILE);

        if (MAX_RECORDS_IN_RAM != null) {
            factory = factory.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
        }

        return ReadUtils.createCommonSAMWriterFromFactory(factory, outputFile, referenceFile, header, preSorted);
    }
}
