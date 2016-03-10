package org.broadinstitute.hellbender.tools.spark;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.codehaus.plexus.util.FileUtils;
import org.seqdoop.hadoop_bam.SplittingBAMIndexer;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

@CommandLineProgramProperties(summary = "Create a HadoopBam splitting index as well as a bam index from a bam file.",
        oneLineSummary = "create a hadoop-bam splitting index" ,
        programGroup = SparkProgramGroup.class)
public final class CreateHadoopBamSplittingIndex extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(CreateHadoopBamSplittingIndex.class);

    public static final String SPLITTING_INDEX_GRANULARITY_LONG_NAME = "splittingIndexGranularity";
    public static final String SPLITTING_INDEX_GRANULARITY_SHORT_NAME = "granularity";
    public static final String CREATE_BAI_LONG_NAME = "createBai";
    public static final String CREATE_BAI_SHORT_NAME = CREATE_BAI_LONG_NAME;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
             doc = "bam file to create a HadoopBam splitting index for",
            optional = false )
    public File inputBam;

    @Argument(fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true,
            optional=true)
    public ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The BAM splitting_index file. If this is unspecified an index will be created with the same name as " +
                    "the input file but with the additional extension " + SplittingBAMIndexer.OUTPUT_FILE_EXTENSION,
            optional = true)
    public File output;

    @Argument(fullName = SPLITTING_INDEX_GRANULARITY_LONG_NAME,
            shortName = SPLITTING_INDEX_GRANULARITY_SHORT_NAME,
            doc = "Splitting index granularity, an entry is created in the index every this many reads.",
            optional = true)
    public int granularity = SplittingBAMIndexer.DEFAULT_GRANULARITY;

    @Argument(fullName = CREATE_BAI_LONG_NAME,
            shortName = CREATE_BAI_SHORT_NAME,
            doc = "set this to create a bai index at the same time as creating a splitting index",
            optional = true)
    public boolean createBai = false;


    @Override
    public Object doWork() {
        if( granularity <= 0) {
            throw new UserException.BadArgumentValue(SPLITTING_INDEX_GRANULARITY_LONG_NAME, Integer.toString(granularity), "Granularity must be > 0");
        }
        final File index = getOutputFile(output, inputBam);
        if(createBai){
            createBaiAndSplittingIndex(inputBam, index, granularity, readValidationStringency);
        } else {
            createOnlySplittingIndex(inputBam, index, granularity);
        }

        return 0;
    }

    private static void createOnlySplittingIndex(final File inputBam, final File index, final int granularity) {
        assertIsBam(inputBam);
        //createBamSplittingIndex(inputBam, getOutputFile(output, inputBam), readValidationStringency, granularity);
        try(BufferedInputStream in = new BufferedInputStream(new FileInputStream(inputBam));
          BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(index))) {
                SplittingBAMIndexer.index(in, out, inputBam.length(), granularity);

        } catch (final IOException e) {
            throw new UserException("Couldn't create splitting index", e);
        }
    }

    private static void createBaiAndSplittingIndex(final File inputBam, final File index, final int granularity, final ValidationStringency readValidationStringency) {
        assertIsBam(inputBam);
        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(readValidationStringency)
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .open(inputBam);
            BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(index))) {
                final SAMFileHeader header = reader.getFileHeader();
                assertBamIsCoordinateSorted(header);
                final SplittingBAMIndexer indexer = new SplittingBAMIndexer(out, granularity);

                final BAMIndexer bamIndexer = new BAMIndexer(IOUtils.replaceExtension(index, BAMIndex.BAMIndexSuffix), header);
                for(final SAMRecord read : reader){
                    indexer.processAlignment(read);
                    bamIndexer.processAlignment(read);
                }
                indexer.finish(inputBam.length());
                bamIndexer.finish();
        } catch (final IOException e) {
            throw new UserException("Couldn't create splitting index", e);
        }
    }

    private static void assertBamIsCoordinateSorted(final SAMFileHeader header) {
        if( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UserException.BadInput("Cannot create a " + BAMIndex.BAMIndexSuffix + " index for a file " +
                    "that isn't coordinate sorted.");
        }
    }


    private static void assertIsBam(final File inputBam) {
        if(!BamFileIoUtils.isBamFile(inputBam)) {
            throw new UserException.BadInput("A splitting index is only relevant for a bam file, but a "
                    + "file with extension "+ FileUtils.getExtension(inputBam.getName()) + " was specified.");
        }
    }

    private static File getOutputFile(final File suggestedOutput, final File input) {
        if(suggestedOutput == null){
            return new File(input.getPath() + SplittingBAMIndexer.OUTPUT_FILE_EXTENSION);
        } else {
            if (!suggestedOutput.getAbsolutePath().endsWith("bam" + SplittingBAMIndexer.OUTPUT_FILE_EXTENSION)){
                logger.warn("Creating a splitting index with an extension that doesn't match "
                        + "bam"+SplittingBAMIndexer.OUTPUT_FILE_EXTENSION + ".  Output file: "+suggestedOutput);
            }
            return suggestedOutput;
        }
    }
}
