package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.*;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Converts a fastq file to an unaligned BAM/SAM format.
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FastQ specification</a> for details.
 * Three fastq versions are supported: FastqSanger, FastqSolexa and FastqIllumina.
 * Input files can be in GZip format (end in .gz).
 */
@CommandLineProgramProperties(
        summary = "Extracts read sequences and qualities from the input fastq file and writes them into the output file in unaligned SAM/BAM format."
                + " Input files can be in GZip format (end in .gz).\n",
        oneLineSummary = "Converts a fastq file to an unaligned SAM/BAM file",
        programGroup = ReadProgramGroup.class
)
public final class FastqToSam extends PicardCommandLineProgram {

    private static final Logger LOG = LogManager.getLogger();

    @Argument(shortName="F1", doc="Input fastq file (optionally gzipped) for single end data, or first read in paired end data.")
    public File FASTQ;

    @Argument(shortName="F2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.", optional=true)
    public File FASTQ2;

    @Argument(shortName="V", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 " +
            "style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled " +
            "scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.", optional = true)
    public FastqQualityFormat QUALITY_FORMAT;

    @Argument(doc="Output SAM/BAM file. ", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName= StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT ;

    @Argument(shortName="RG", doc="Read group name")
    public String READ_GROUP_NAME = "A";

    @Argument(shortName="SM", doc="Sample name to insert into the read group header")
    public String SAMPLE_NAME;

    @Argument(shortName="LB", doc="The library name to place into the LB attribute in the read group header", optional=true)
    public String LIBRARY_NAME;

    @Argument(shortName="PU", doc="The platform unit (often run_barcode.lane) to insert into the read group header", optional=true)
    public String PLATFORM_UNIT;

    @Argument(shortName="PL", doc="The platform type (e.g. illumina, solid) to insert into the read group header", optional=true)
    public String PLATFORM;

    @Argument(shortName="CN", doc="The sequencing center from which the data originated", optional=true)
    public String SEQUENCING_CENTER;

    @Argument(shortName = "PI", doc = "Predicted median insert size, to insert into the read group header", optional = true)
    public Integer PREDICTED_INSERT_SIZE;

    @Argument(shortName = "PG", doc = "Program group to insert into the read group header.", optional = true)
    public String PROGRAM_GROUP;

    @Argument(shortName = "PM", doc = "Platform model to insert into the group header " +
            "(free-form text providing further details of the platform/technology used)", optional = true)
    public String PLATFORM_MODEL;

    @Argument(doc="Comment(s) to include in the merged output file's header.", optional=true, shortName="CO")
    public List<String> COMMENT = new ArrayList<>();

    @Argument(shortName = "DS", doc = "Inserted into the read group header", optional = true)
    public String DESCRIPTION;

    @Argument(shortName = "DT", doc = "Date the run was produced, to insert into the read group header", optional = true)
    public Iso8601Date RUN_DATE;

    @Argument(shortName="SO", doc="The sort order for the output sam/bam file.")
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.queryname;

    @Argument(doc="Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.")
    public int MIN_Q = 0;

    @Argument(doc="Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.")
    public int MAX_Q = SAMUtils.MAX_PHRED_SCORE;

    @Argument(doc="If true and this is an unpaired fastq any occurance of '/1' will be removed from the end of a read name.")
    public Boolean STRIP_UNPAIRED_MATE_NUMBER = false;

    @Argument(doc="Allow (and ignore) empty lines")
    public Boolean ALLOW_AND_IGNORE_EMPTY_LINES = false;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    /**
     * Looks at fastq input(s) and attempts to determine the proper quality format
     *
     * Closes the reader(s) by side effect
     *
     * @param reader1 The first fastq input
     * @param reader2 The second fastq input, if necessary. To not use this input, set it to null
     * @param expectedQuality If provided, will be used for sanity checking. If left null, autodetection will occur
     */
    public static FastqQualityFormat determineQualityFormat(final FastqReader reader1, final FastqReader reader2, final FastqQualityFormat expectedQuality) {
        final QualityEncodingDetector detector = new QualityEncodingDetector();

        if (reader2 == null) {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1);
        } else {
            detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, reader1, reader2);
            reader2.close();
        }

        reader1.close();

        final FastqQualityFormat qualityFormat =  detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, expectedQuality);
        if (detector.isDeterminationAmbiguous()) {
            LOG.warn("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.");
        }
        LOG.info(String.format("Auto-detected quality format as: %s.", qualityFormat));

        return qualityFormat;
    }

    /* Simply invokes the right method for unpaired or paired data. */
    protected Object doWork() {
        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(OUTPUT);

        FastqReader reader = fileToFastqReader(FASTQ);

        FastqReader reader2 = null;
        if (FASTQ2 != null) {
            IOUtil.assertFileIsReadable(FASTQ2);
            reader2 = fileToFastqReader(FASTQ2);
        }

        QUALITY_FORMAT = FastqToSam.determineQualityFormat(reader, reader2, QUALITY_FORMAT);

        final SAMFileHeader header = createSamFileHeader();
        try (final SAMFileWriter writer = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, header, false)) {

            reader = fileToFastqReader(FASTQ);
            if (FASTQ2 != null) {
                reader2 = fileToFastqReader(FASTQ2);
            }
            makeItSo(reader, reader2, writer);

            reader.close();

            if (reader2 != null) {
                reader2.close();
            }
        }
        return null;
    }

    /**
     * Handles the FastqToSam execution on the FastqReader(s).
     *
     * In some circumstances it might be useful to circumvent the command line based instantiation of this
     * class, however note that there is no handholding or guardrails to running in this manner.
     *
     * It is the caller's responsibility to close the reader(s)
     *
     * @param reader1 The FastqReader for the first fastq file
     * @param reader2 The second FastqReader if applicable. Pass in null if only using a single reader
     * @param writer The SAMFileWriter where the new SAM file is written
     *
     */
    public void makeItSo(final FastqReader reader1, final FastqReader reader2, final SAMFileWriter writer) {
        final int readCount = (reader2 == null) ?  doUnpaired(reader1, writer) : doPaired(reader1, reader2, writer);
        LOG.info("Processed " + readCount + " fastq reads");
    }

    /** Creates a simple SAM file from a single fastq file. */
    protected int doUnpaired(final FastqReader freader, final SAMFileWriter writer) {
        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader.hasNext()  ; readCount++) {
            final FastqRecord frec = freader.next();
            final SAMRecord srec = createSamRecord(writer.getFileHeader(), getReadName(frec.getReadHeader(), false) , frec, false) ;
            srec.setReadPairedFlag(false);
            writer.addAlignment(srec);
            progress.record(srec);
        }

        writer.close();
        return readCount;
    }

    /** More complicated method that takes two fastq files and builds pairing information in the SAM. */
    protected int doPaired(final FastqReader freader1, final FastqReader freader2, final SAMFileWriter writer) {
        int readCount = 0;
        final ProgressLogger progress = new ProgressLogger(LOG);
        for ( ; freader1.hasNext() && freader2.hasNext() ; readCount++) {
            final FastqRecord frec1 = freader1.next();
            final FastqRecord frec2 = freader2.next();

            final String frec1Name = getReadName(frec1.getReadHeader(), true);
            final String frec2Name = getReadName(frec2.getReadHeader(), true);
            final String baseName = getBaseName(frec1Name, frec2Name, freader1, freader2);

            final SAMRecord srec1 = createSamRecord(writer.getFileHeader(), baseName, frec1, true) ;
            srec1.setFirstOfPairFlag(true);
            srec1.setSecondOfPairFlag(false);
            writer.addAlignment(srec1);
            progress.record(srec1);

            final SAMRecord srec2 = createSamRecord(writer.getFileHeader(), baseName, frec2, true) ;
            srec2.setFirstOfPairFlag(false);
            srec2.setSecondOfPairFlag(true);
            writer.addAlignment(srec2);
            progress.record(srec2);
        }

        writer.close();

        if (freader1.hasNext() || freader2.hasNext()) {
            throw new UserException("Input paired fastq files must be the same length");
        }

        return readCount;
    }

    private FastqReader fileToFastqReader(final File file) {
        return new FastqReader(file, ALLOW_AND_IGNORE_EMPTY_LINES);
    }

    private SAMRecord createSamRecord(final SAMFileHeader header, final String baseName, final FastqRecord frec, final boolean paired) {
        final SAMRecord srec = new SAMRecord(header);
        srec.setReadName(baseName);
        srec.setReadString(frec.getReadString());
        srec.setReadUnmappedFlag(true);
        srec.setAttribute(ReservedTagConstants.READ_GROUP_ID, READ_GROUP_NAME);
        final byte[] quals = StringUtil.stringToBytes(frec.getBaseQualityString());
        convertQuality(quals, QUALITY_FORMAT);
        for (final byte qual : quals) {
            final int uQual = qual & 0xff;
            if (uQual < MIN_Q || uQual > MAX_Q) {
                throw new GATKException("Base quality " + uQual + " is not in the range " + MIN_Q + ".." +
                        MAX_Q + " for read " + frec.getReadHeader());
            }
        }
        srec.setBaseQualities(quals);

        if (paired) {
            srec.setReadPairedFlag(true);
            srec.setMateUnmappedFlag(true);
        }
        return srec ;
    }

    /** Creates a simple header with the values provided on the command line. */
    public SAMFileHeader createSamFileHeader() {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(this.READ_GROUP_NAME);
        rgroup.setSample(this.SAMPLE_NAME);
        if (this.LIBRARY_NAME != null) rgroup.setLibrary(this.LIBRARY_NAME);
        if (this.PLATFORM != null) rgroup.setPlatform(this.PLATFORM);
        if (this.PLATFORM_UNIT != null) rgroup.setPlatformUnit(this.PLATFORM_UNIT);
        if (this.SEQUENCING_CENTER != null) rgroup.setSequencingCenter(SEQUENCING_CENTER);
        if (this.PREDICTED_INSERT_SIZE != null) rgroup.setPredictedMedianInsertSize(PREDICTED_INSERT_SIZE);
        if (this.DESCRIPTION != null) rgroup.setDescription(this.DESCRIPTION);
        if (this.RUN_DATE != null) rgroup.setRunDate(this.RUN_DATE);
        if (this.PLATFORM_MODEL != null) rgroup.setPlatformModel(this.PLATFORM_MODEL);
        if (this.PROGRAM_GROUP != null) rgroup.setProgramGroup(this.PROGRAM_GROUP);

        final SAMFileHeader header = new SAMFileHeader();
        header.addReadGroup(rgroup);

        for (final String comment : COMMENT) {
            header.addComment(comment);
        }

        header.setSortOrder(this.SORT_ORDER);
        return header ;
    }

    /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
    void convertQuality(final byte[] quals, final FastqQualityFormat version) {
        switch (version)  {
            case Standard:
                SAMUtils.fastqToPhred(quals);
                break ;
            case Solexa:
                solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
                break ;
            case Illumina:
                solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
                break ;
        }
    }

    /** Returns read baseName and asserts correct pair read name format:
     * <ul>
     * <li> Paired reads must either have the exact same read names or they must contain at least one "/"
     * <li> and the First pair read name must end with "/1" and second pair read name ends with "/2"
     * <li> The baseName (read name part before the /) must be the same for both read names
     * <li> If the read names are exactly the same but end in "/2" or "/1" then an exception will be thrown
     * </ul>
     */
    String getBaseName(final String readName1, final String readName2, final FastqReader freader1, final FastqReader freader2) {
        String [] toks = getReadNameTokens(readName1, 1, freader1);
        final String baseName1 = toks[0] ;
        final String num1 = toks[1] ;

        toks = getReadNameTokens(readName2, 2, freader2);
        final String baseName2 = toks[0] ;
        final String num2 = toks[1];

        if (!baseName1.equals(baseName2)) {
            throw new UserException(String.format("In paired mode, read name 1 (%s) does not match read name 2 (%s)", baseName1,baseName2));
        }

        final boolean num1Blank = StringUtil.isBlank(num1);
        final boolean num2Blank = StringUtil.isBlank(num2);
        if (num1Blank || num2Blank) {
            if(!num1Blank) throw new UserException(error(freader1,"Pair 1 number is missing (" +readName1+ "). Both pair numbers must be present or neither."));       //num1 != blank and num2   == blank
            else if(!num2Blank) throw new UserException(error(freader2, "Pair 2 number is missing (" +readName2+ "). Both pair numbers must be present or neither.")); //num1 == blank and num =2 != blank
        } else {
            if (!num1.equals("1")) throw new UserException(error(freader1,"Pair 1 number must be 1 ("+readName1+")"));
            if (!num2.equals("2")) throw new UserException(error(freader2,"Pair 2 number must be 2 ("+readName2+")"));
        }

        return baseName1 ;
    }

    /** Breaks up read name into baseName and number separated by the last / */
    private String [] getReadNameTokens(final String readName, final int pairNum, final FastqReader freader) {
        if(readName.equals("")) throw new UserException(error(freader,"Pair read name "+pairNum+" cannot be empty: "+readName));

        final int idx = readName.lastIndexOf("/");
        final String result[] = new String[2];

        if (idx == -1) {
            result[0] = readName;
            result[1] = null;
        } else {
            result[1] = readName.substring(idx+1, readName.length()); // should be a 1 or 2

            if(!result[1].equals("1") && !result[1].equals("2")) {    //if not a 1 or 2 then names must be identical
                result[0] = readName;
                result[1] = null;
            }
            else {
                result[0] = readName.substring(0,idx); // baseName
            }
        }

        return result ;
    }

    /** Little utility to give error messages corresponding to line numbers in the input files. */
    private String error(final FastqReader freader, final String str) {
        return str +" at line "+freader.getLineNumber() +" in file "+freader.getFile().getAbsolutePath();
    }

    // Read names cannot contain blanks
    private String getReadName(final String fastqHeader, final boolean paired) {
        final int idx = fastqHeader.indexOf(" ");
        String readName = (idx == -1) ? fastqHeader : fastqHeader.substring(0,idx);

        // NOTE: the while loop isn't necessarily the most efficient way to handle this but we don't
        // expect this to ever happen more than once, just trapping pathological cases
        while (STRIP_UNPAIRED_MATE_NUMBER && !paired && (readName.endsWith("/1") || readName.endsWith("/2"))) {
            // If this is an unpaired run we want to make sure that "/1" isn't tacked on the end of the read name,
            // as this can cause problems down the road in MergeBamAlignment
            readName = readName.substring(0, readName.length() - 2);
        }

        return readName;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (MIN_Q < 0) return new String[]{"MIN_Q must be >= 0"};
        if (MAX_Q > SAMUtils.MAX_PHRED_SCORE) return new String[]{"MAX_Q must be <= " + SAMUtils.MAX_PHRED_SCORE};
        return null;
    }
}
