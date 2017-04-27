package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.File;
import java.io.IOException;
import java.nio.channels.SeekableByteChannel;
import java.util.*;
import java.util.function.Function;

/**
 * This tool imports GVCFs to GenomicsDB. To run this tool,
 * 1. A single interval must be provided
 * 2. The tool accepts multiple GVCFs each of which must contain data
 *    for one sample
 * 3. The path to the GenomicsDB workspace must be specified
 * 4. User may optionally specify paths to which to write JSON files
 *
 * To read data from GenomicsDB, use the query interface GenomicsDBFeatureReader
 */
@CommandLineProgramProperties(
    summary = "Import VCFs to GenomicsDB",
    oneLineSummary = "Import VCFs to GenomicsDB",
    programGroup = VariantProgramGroup.class
)
public final class GenomicsDBImport extends GATKTool {

    private static final long DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE = 16*1024L;
    private static final long DEFAULT_SEGMENT_SIZE = 1048576L;

    public static final String WORKSPACE_ARG_NAME = "genomicsDBWorkspace";
    private static final String SEGMENT_SIZE_ARG_NAME = "genomicsDBSegmentSize";
    private static final String OVERWRITE_WORKSPACE_NAME = "overwriteExistingGenomicsDBWorkspace";

    private static final String VCF_BUFFER_SIZE_ARG_NAME = "genomicsDBVCFBufferSize";
    private static final String ARRAY_ARG_NAME = "genomicsDBArray";
    private static final String VID_MAP_FILE_ARG_NAME = "genomicsDBVidMapFile";

    private static final String CALLSET_MAP_FILE_ARG_NAME = "genomicsDBCallsetMapFile";

    @Argument(fullName = WORKSPACE_ARG_NAME,
              shortName = WORKSPACE_ARG_NAME,
              doc = "Workspace for GenomicsDB. Has to be a POSIX file system path")
    private String workspace;

//    @Argument(fullName = ARRAY_ARG_NAME,
//              shortName = ARRAY_ARG_NAME,
//              doc = "TileDB array name used by GenomicsDB. Defaults to " + ARRAY_ARG_NAME,
//              optional = true)
    private final String arrayName = GenomicsDBConstants.DEFAULT_ARRAY_NAME;

    @Argument(fullName = SEGMENT_SIZE_ARG_NAME,
              shortName = SEGMENT_SIZE_ARG_NAME,
              doc = "Buffer size in bytes allocated for GenomicsDB attributes during " +
                    "import. Should be large enough to hold data from one site. " +
                    " Defaults to " + DEFAULT_SEGMENT_SIZE,
              optional = true)
    private long segmentSize = DEFAULT_SEGMENT_SIZE;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              doc = "GVCF files to be imported to GenomicsDB. Each file must contain" +
                    "data for only a single sample")
    private List<String> variantFilePaths;

    @Argument(fullName = VCF_BUFFER_SIZE_ARG_NAME,
              shortName = VCF_BUFFER_SIZE_ARG_NAME,
              doc = "Buffer size in bytes to store variant contexts." +
                    " Larger values are better as smaller values cause frequent disk writes." +
                    " Defaults to " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE,
              optional = true)
    private long vcfBufferSizePerSample = DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE;

//    @Argument(fullName = CALLSET_MAP_FILE_ARG_NAME,
//              shortName = CALLSET_MAP_FILE_ARG_NAME,
//              doc = "Path to callset JSON file to be created. " +
//                    "The file contains row mappings to sample names."
//                    + " Defaults to " + DEFAULT_CALLSETMAP_FILE_NAME +
//                    " file is written to GenomicsDB workspace",
//              optional = true)
    private final String callsetMapJSONFileName = GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME;

//    @Argument(fullName = VID_MAP_FILE_ARG_NAME,
//              shortName = VID_MAP_FILE_ARG_NAME,
//              doc = "Path to vid map JSON file to be created. Includes contig maps and INFO/FORMAT/FILTER" +
//                    "fields from headers. Defaults to " + DEFAULT_VIDMAP_FILE_NAME +
//                    " file is written to GenomicsDB workspace",
//              optional = true)
    private final String vidMapJSONFileName = GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME;

    @Argument(fullName = OVERWRITE_WORKSPACE_NAME,
              shortName = OVERWRITE_WORKSPACE_NAME,
              doc = "Will overwrite given workspace if it exists. " +
                    "Otherwise a new workspace is created. " +
                    "Defaults to false",
              optional = true)
    private Boolean overwriteExistingWorkspace = false;

    @Override
    public boolean requiresIntervals() { return true; }

    @Override
    public boolean requiresReads() { return false; }

    @Override
    public boolean requiresReference() {
      return false;
    }

    // The GenomicsDB importer object
    private GenomicsDBImporter importer = null;

    // Intervals from command line (singleton for now)
    private List<ChromosomeInterval> intervals;

    private final Map<String, FeatureReader<VariantContext>> sampleToVCFMap = new HashMap<>();

    // Needed as smartMergeHeaders() returns a set of VCF header lines
    private Set<VCFHeaderLine> mergedHeader = null;

    // VCFHeader object created from mergedHeader
    private VCFHeader mergedVCFHeader = null;

    /**
     * Before traversal starts, create the feature readers
     * for all the input GVCFs, create the merged header and
     * initialize the interval
     */
    @Override
    public void onStartup() {

        for (final String variantPath : variantFilePaths) {

            final String variantURI = IOUtils.getPath(variantPath).toAbsolutePath().toUri().toString();
            final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = (cloudPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudPrefetchBuffer, is) : Function.identity());
            final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = (cloudIndexPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudIndexPrefetchBuffer, is) : Function.identity());
            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), true, cloudWrapper, cloudIndexWrapper);

            // A GVCF file must contain only one sample, throw an exception otherwise
            if (((VCFHeader)reader.getHeader()).getNGenotypeSamples() != 1) {
                throw new UserException("Input GVCF: " + variantPath + " should contain data for one sample");
            }

            final String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
            sampleToVCFMap.put(sampleName, reader);
        }

        createMergedHeader();

        initializeIntervals();

        super.onStartup();
    }

    /**
     * Before traversal, fix configuration parameters and initialize
     * GenomicsDB. Hard-coded to handle only VCF files and headers
     */
    @Override
    public void onTraversalStart() {

        File workspaceDir = overwriteOrCreateWorkspace();

        // If provided buffer size is less 1KB it is considered to be too small
        // Either use default 16KB or any value larger than 10KB
        if (vcfBufferSizePerSample < 1024L) {
            logger.warn("Buffer size per column partition per sample is too small." +
                " Either use default value " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE +
                " or larger values than 10KB");
        }

        final File vidMapJSONFile = (vidMapJSONFileName.equals(GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME)) ?
            new File(workspaceDir + "/" + vidMapJSONFileName) :
            new File(vidMapJSONFileName);

        final File callsetMapJSONFile = (callsetMapJSONFileName.equals(GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME)) ?
            new File (workspaceDir + "/" + callsetMapJSONFileName) :
            new File(callsetMapJSONFileName);

        logger.info("Vid Map JSON file will be written to " + vidMapJSONFile);
        logger.info("Callset Map JSON file will be written to " + callsetMapJSONFile);
        
        final long variantContextBufferSize = vcfBufferSizePerSample * sampleToVCFMap.size();

        logger.info("Writing data to array - " + workspace + "/" + arrayName);

        try {
            logger.info("GenomicsDB import batch started");
            importer = new GenomicsDBImporter(sampleToVCFMap, mergedHeader, intervals.get(0), workspace,
                                              arrayName, variantContextBufferSize, segmentSize,
                                              vidMapJSONFile.getAbsolutePath(),
                                              callsetMapJSONFile.getAbsolutePath());
        } catch (IOException e) {
            throw new UserException("Error initializing GenomicsDBImporter", e);
        }
    }

    /**
     * A complete traversal from start to finish. This method will import all samples
     * specified in the input GVCF files.
     */
    @Override
    public void traverse() {
        try {
            importer.importBatch();
        } catch (IOException e) {
            throw new UserException("GenomicsDB importBatch raised exception", e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Import completed!");
        return importer.isDone();
    }

    @Override
    public void closeTool() {
        for (Map.Entry<String, FeatureReader<VariantContext>> reader : sampleToVCFMap.entrySet()) {
            try {
               reader.getValue().close();
            } catch (IOException e) {
                throw new GATKException("FeatureReader close() failed for " + reader.getKey(), e);
            }
        }
    }

    /**
     * Input argument "overwriteExistingWorkspace" defaults to false.
     * The tool creates a new workspace if it doesn't exist. Deletes
     * an existing workspace if argument is true
     *
     * @return  The workspace directory
     */
    private File overwriteOrCreateWorkspace() {
        File workspaceDir = new File(workspace);

        if (overwriteExistingWorkspace) {
            IOUtils.tryDelete(workspaceDir);
        }

        if (!workspaceDir.exists()) {
            int ret = GenomicsDBImporter.createTileDBWorkspace(workspaceDir.getAbsolutePath());

            if (ret > 0) {
                checkIfValidWorkspace(workspaceDir);
                logger.info("Importing data to GenomicsDB workspace: " + workspaceDir);
            } else if (ret < 0) {
                throw new UserException("Error creating GenomicsDB workspace: " + workspaceDir);
            }
        } else {
            // Check whether its a valid workspace
            checkIfValidWorkspace(workspaceDir);
        }
        return workspaceDir;
    }

    private void checkIfValidWorkspace(File workspaceDir) {
        File tempFile = new File(workspaceDir.getAbsolutePath() + "/__tiledb_workspace.tdb");
        if (!tempFile.exists()) {
            throw new UserException(workspaceDir.getAbsolutePath() + " is not a valid GenomicsDB workspace");
        }
    }

    /**
     * Create merged header from input GVCFs
     */
    private void createMergedHeader() {
        List<VCFHeader> headers = new ArrayList<>();
        for (Map.Entry<String, FeatureReader<VariantContext>> sampleToVCF : sampleToVCFMap.entrySet()) {
            FeatureReader<VariantContext> reader = sampleToVCF.getValue();
            VCFHeader vcfHeader = (VCFHeader)reader.getHeader();
            headers.add((VCFHeader)reader.getHeader());
        }

        mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
        mergedVCFHeader = new VCFHeader(mergedHeader);
    }

    /**
     * Loads our intervals using the best available sequence
     * dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
     * to parse/verify them. Does nothing if no intervals were specified.
     */
    private void initializeIntervals() {
        if (intervalArgumentCollection.intervalsSpecified()) {

            final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();
            if (intervalDictionary == null) {
                throw new UserException("We require at least one input source that " +
                    "has a sequence dictionary (reference or reads) when intervals are specified");
            }

            intervals = new ArrayList<>();

            List<SimpleInterval> simpleIntervalList =
                intervalArgumentCollection.getIntervals(intervalDictionary);

            if (simpleIntervalList.size() > 1) {
                throw new UserException("More than one interval specified. The tool takes only one");
            }

            for (SimpleInterval simpleInterval : simpleIntervalList) {
                intervals.add(new ChromosomeInterval(simpleInterval.getContig(),
                  simpleInterval.getStart(), simpleInterval.getEnd()));
            }
        } else {
            throw new UserException("No intervals specified");
        }
    }

    /**
     * Overriding getBestAvailableSequenceDictionary() to prefer the mergedVCFHeader's
     * sequence directory, if present, over any other dictionaries
     *
     * @return  Sequence directory from merged header, or super.getBestAvailableSequenceDictionary()
     *          if none
     */
    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        final SAMSequenceDictionary sequenceDictionary = mergedVCFHeader.getSequenceDictionary();
        if (sequenceDictionary == null) {
            return super.getBestAvailableSequenceDictionary();
        } else {
            return sequenceDictionary;
        }
    }
}

