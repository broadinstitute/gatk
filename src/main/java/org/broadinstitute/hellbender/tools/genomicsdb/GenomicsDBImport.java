package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBCallsetsMapProto;
import com.intel.genomicsdb.GenomicsDBImportConfiguration;
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
import java.io.FileNotFoundException;
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
    private static final int DEFAULT_ZERO_BATCH_SIZE = 0;

    public static final String WORKSPACE_ARG_NAME = "genomicsDBWorkspace";
    private static final String SEGMENT_SIZE_ARG_NAME = "genomicsDBSegmentSize";
    private static final String OVERWRITE_WORKSPACE_NAME = "overwriteExistingGenomicsDBWorkspace";

    private static final String VCF_BUFFER_SIZE_ARG_NAME = "genomicsDBVCFBufferSize";
    private static final String ARRAY_ARG_NAME = "genomicsDBArray";
    private static final String VID_MAP_FILE_ARG_NAME = "genomicsDBVidMapFile";

    private static final String CALLSET_MAP_FILE_ARG_NAME = "genomicsDBCallsetMapFile";
    private static final String BATCHSIZE_ARG_NAME = "batchSize";
    private static final String CONSOLIDATE_ARG_NAME = "consolidate";

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
    private List<String> variantPaths;

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

    @Argument(fullName = BATCHSIZE_ARG_NAME,
              shortName = BATCHSIZE_ARG_NAME,
              doc = "Batch size controls the number of samples for which readers are open at once " +
                    "and therefore provides a way to minimize memory consumption. However, it can take longer to complete. " +
                    "Use the consolidate flag if more than a hundred batches were used. This will improve feature read time. " +
                    "batchSize=0 means no batching (i.e. readers for all samples will be opened at once) " +
                    "Defaults to " + DEFAULT_ZERO_BATCH_SIZE,
              optional = true)
    private int batchSize = DEFAULT_ZERO_BATCH_SIZE;

    @Argument(fullName = CONSOLIDATE_ARG_NAME,
              shortName = CONSOLIDATE_ARG_NAME,
              doc = "Boolean flag to enable consolidation. If importing data in batches, a new fragment is created for " +
                    "each batch. In case thousands of fragments are created, GenomicsDB feature readers will try " +
                    "to open ~20x as many files. Also, internally GenomicsDB would consume more memory to maintain " +
                    "bookkeeping data from all fragments. Use this flag to merge all fragments into one. " +
                    "Merging can potentially improve read performance, however overall benefit might not be noticeable " +
                    "as the top Java layers has significantly higher overheads. This flag have no effect if only one " +
                    "batch is used. Defaults to false",
              optional = true)
    private Boolean doConsolidation = false;

    @Override
    public boolean requiresIntervals() { return true; }

    @Override
    public boolean requiresReads() { return false; }

    @Override
    public boolean requiresReference() {
      return false;
    }

    // Intervals from command line (singleton for now)
    private List<ChromosomeInterval> intervals;

    // List of all sample names
    private List<String> sampleNames = new ArrayList<>();

    // Linked hash map between sample names and corresponding GVCF file name
    private Map<String, String> sampleNameToVcfUri = new LinkedHashMap<>();

    // Needed as smartMergeHeaders() returns a set of VCF header lines
    private Set<VCFHeaderLine> mergedHeader = null;

    // VCFHeader created from the merged header
    private VCFHeader mergedVCFHeader;

    // Path to vidmap file to be written by GenomicsDBImporter
    private File vidMapJSONFile;

    // Path to callsetmap file to be written by GenomicsDBImporter
    private File callsetMapJSONFile;

    // GenomicsDB callset map protobuf structure containing all callset names
    // used to write the callset json file on traversal success
    private GenomicsDBCallsetsMapProto.CallsetMappingPB callsetMappingPB;

    /**
     * Before traversal starts, create the feature readers
     * for all the input GVCFs, create the merged header and
     * initialize the interval
     */
    @Override
    public void onStartup() {

        List<VCFHeader> headers = new ArrayList<>(variantPaths.size());

        for (String variantPath : variantPaths) {

            final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromVCFUri(variantPath);
            VCFHeader header = (VCFHeader)reader.getHeader();

            // A GVCF file must contain only one sample, throw an exception otherwise
            if (header.getNGenotypeSamples() != 1) {
                throw new UserException("Input GVCF: " + variantPath + " should contain data for one sample");
            }
            String sampleName = header.getGenotypeSamples().get(0);
            headers.add(header);
            sampleNames.add(sampleName);
            if (sampleNameToVcfUri.put(sampleName, variantPath) != null) {
                throw new UserException("Duplicate sample " + sampleName);
            }

            // Import occurs in batches. By closing all readers here we release resources
            // and preserve memory for later stages where readers are allocated
            // for samples in a given batch
            try {
                reader.close();
            } catch (IOException e) {
                throw new GATKException("FeatureReader close() failed for " + sampleName, e);
            }
        }

        mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
        mergedVCFHeader = new VCFHeader(mergedHeader);

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
            throw new UserException("Buffer size per column partition per sample is too small." +
                    " Either use default value " + DEFAULT_VCF_BUFFER_SIZE_PER_SAMPLE +
                    " or larger values than 10KB");
        }

        vidMapJSONFile = (vidMapJSONFileName.equals(GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME)) ?
            new File(workspaceDir + "/" + vidMapJSONFileName) :
            new File(vidMapJSONFileName);

        callsetMapJSONFile = (callsetMapJSONFileName.equals(GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME)) ?
            new File (workspaceDir + "/" + callsetMapJSONFileName) :
            new File(callsetMapJSONFileName);

        logger.info("Vid Map JSON file will be written to " + vidMapJSONFile);
        logger.info("Callset Map JSON file will be written to " + callsetMapJSONFile);
        logger.info("Importing to array - " + workspace + "/" + arrayName);

        // Passing in false here so that sample names will be sorted.
        // This is needed for consistent ordering across partitions/machines
        callsetMappingPB = GenomicsDBImporter.generateSortedCallSetMap(sampleNames, false);
    }

    /**
     * A complete traversal from start to finish. This method will import all samples
     * specified in the input GVCF files.
     */
    @Override
    public void traverse() {

        int updatedBatchSize = (batchSize == DEFAULT_ZERO_BATCH_SIZE) ? sampleNames.size() : batchSize;
        int sampleCount = sampleNames.size();
        int totalBatchCount = (sampleCount/updatedBatchSize) + (sampleCount%updatedBatchSize==0 ? 0 : 1);

        GenomicsDBImporter importer;

        for (int i = 0, batchCount = 1; i < sampleCount; i += updatedBatchSize, ++batchCount) {

            Map<String, FeatureReader<VariantContext>> sampleToReaderMap =
                getFeatureReaders(sampleNames, sampleNameToVcfUri, updatedBatchSize, sampleCount, i);

            logger.info("Importing batch " + batchCount + " with " + sampleToReaderMap.size() + " samples");
            final long variantContextBufferSize = vcfBufferSizePerSample * sampleToReaderMap.size();

            GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
                    createImportConfiguration(workspace, arrayName,
                            variantContextBufferSize, segmentSize,
                            i, (i+updatedBatchSize-1));

            try {
                importer = new GenomicsDBImporter(sampleToReaderMap, mergedHeader, intervals.get(0), importConfiguration);
            } catch (IOException e) {
                throw new UserException("Error initializing GenomicsDBImporter in batch " + batchCount, e);
            }

            try {
                importer.importBatch();
            } catch (IOException e) {
                throw new UserException("GenomicsDB import failed in batch " + batchCount, e);
            }

            closeReaders(sampleToReaderMap);
            logger.info("Done importing batch " + batchCount + "/" + totalBatchCount);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (batchSize==DEFAULT_ZERO_BATCH_SIZE) {
            logger.info("Import completed!");
        } else {
            logger.info("Import of all batches to GenomicsDB completed!");
        }

        // Write the vid and callset map JSON files
        try {
            GenomicsDBImporter.writeVidMapJSONFile(vidMapJSONFile.getAbsolutePath(), mergedHeader);
        } catch (FileNotFoundException fe) {
            throw new UserException("Unable to write vid map JSON file " + vidMapJSONFile.getAbsolutePath(), fe);
        }
        try {
            GenomicsDBImporter.writeCallsetMapJSONFile(callsetMapJSONFile.getAbsolutePath(), callsetMappingPB);
        } catch (FileNotFoundException fe) {
            throw new UserException("Unable to write callset map JSON file " + callsetMapJSONFile.getAbsolutePath(), fe);
        }

        if (doConsolidation) {
            logger.info("GenomicsDB consolidation started");
            GenomicsDBImporter.consolidateTileDBArray(workspace, arrayName);
            logger.info("GenomicsDB consolidation completed");
        }

        return true;
    }

    /**
     * Method to create feature readers for input files or GCS URLs
     * in the current batch
     *
     * @param sampleNames  List of all sample names
     * @param sampleNameToFileName  Sample name to file name mapping
     * @param batchSize  Current batch size
     * @param sampleCount  total samples in this import
     * @param lowerSampleIndex  0-based Lower bound of sample index -- inclusive
     * @return  Feature readers to be imported in the current batch
     */
    private Map<String, FeatureReader<VariantContext>> getFeatureReaders(List<String> sampleNames, Map<String, String> sampleNameToFileName, int batchSize, int sampleCount, int lowerSampleIndex) {
        Map<String, FeatureReader<VariantContext>> sampleToReaderMap = new LinkedHashMap<>();

        for(int i = lowerSampleIndex; i < sampleCount && i < lowerSampleIndex+batchSize; ++i) {
            final String sampleName = sampleNames.get(i);
            assert sampleNameToFileName.containsKey(sampleName);

            final AbstractFeatureReader<VariantContext, LineIterator> reader = getReaderFromVCFUri(sampleNameToFileName.get(sampleName));

            assert sampleName.equals(((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0));
            sampleToReaderMap.put(sampleName, reader);
        }
        return sampleToReaderMap;
    }

    /**
     * Creates a feature reader object from a given VCF URI (can also be
     * a local file path) and returns it
     * @param variantPath  URI or file path
     * @return  Feature reader
     */
    private AbstractFeatureReader<VariantContext, LineIterator> getReaderFromVCFUri(String variantPath) {
        final String variantURI = IOUtils.getPath(variantPath).toAbsolutePath().toUri().toString();
        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = (cloudPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudPrefetchBuffer, is) : Function.identity());
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = (cloudIndexPrefetchBuffer > 0 ? is -> SeekableByteChannelPrefetcher.addPrefetcher(cloudIndexPrefetchBuffer, is) : Function.identity());
        return AbstractFeatureReader.getFeatureReader(variantURI, null, new VCFCodec(), true, cloudWrapper, cloudIndexWrapper);
    }

    /**
     * Creates a GenomicsDB configuration data structure
     * instead of sending a long list of parameters to the constructor call
     *
     * @param workspace  GenomicsDB workspace
     * @param arrayName  GenomicsDB array
     * @param variantContextBufferSize  Buffer size to store VCF records for all samples
     * @param segmentSize  Buffer size to store columnar data to be serialized to disk
     * @param lbSampleIndex  Lower bound of sample index -- inclusive (0-based)
     * @param ubSampleIndex  Upper bound of sample index -- inclusive (0-based)
     * @return  GenomicsDB import configuration object
     */
    private GenomicsDBImportConfiguration.ImportConfiguration createImportConfiguration(
        final String workspace,
        final String arrayName,
        final long variantContextBufferSize,
        final long segmentSize,
        final long lbSampleIndex,
        final long ubSampleIndex) {

        GenomicsDBImportConfiguration.Partition.Builder pBuilder =
            GenomicsDBImportConfiguration.Partition.newBuilder();

        // Since, there is one partition for this import, the
        // begin column partition index is 0
        GenomicsDBImportConfiguration.Partition partition =
            pBuilder
                .setWorkspace(workspace)
                .setArray(arrayName)
                .setBegin(0)
                .build();

        GenomicsDBImportConfiguration.GATK4Integration.Builder gBuilder =
            GenomicsDBImportConfiguration.GATK4Integration.newBuilder();

        GenomicsDBImportConfiguration.GATK4Integration gatk4Parameters =
            gBuilder
                .setLowerSampleIndex(lbSampleIndex)
                .setUpperSampleIndex(ubSampleIndex)
                .build();

        GenomicsDBImportConfiguration.ImportConfiguration.Builder cBuilder =
            GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();

        return cBuilder
                .addColumnPartitions(0, partition)
                .setGatk4IntegrationParameters(gatk4Parameters)
                .setSizePerColumnPartition(variantContextBufferSize)
                .setSegmentSize(segmentSize)
                .build();
    }

    /**
     * Close all readers in the current batch
     *
     * @param sampleToReaderMap  Map of sample names to readers
     */
    private void closeReaders(Map<String, FeatureReader<VariantContext>> sampleToReaderMap) {
        for (Map.Entry<String, FeatureReader<VariantContext>> reader : sampleToReaderMap.entrySet()) {
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

