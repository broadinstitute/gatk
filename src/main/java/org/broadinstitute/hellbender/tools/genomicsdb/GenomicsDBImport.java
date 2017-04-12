package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBImporter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.VariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.GenomicsDBProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;

/**
 * This tool imports GVCFs to GenomicsDB.
 * User must specify a loader JSON configuration file,
 * callsets JSON file containing the stream names for samples,
 * and a stream JSON file containing the filenames of the streams
 *
 * Use the GenomicsDB query interface for reading variants from
 * the database
 */
@CommandLineProgramProperties(
  summary = "Import VCFs to GenomicsDB",
  oneLineSummary = "Import VCFs to GenomicsDB",
  programGroup = GenomicsDBProgramGroup.class
)
public final class GenomicsDBImport extends GATKTool {

  private static final long serialVersionUID = 1L;
  private final long DEFAULT_SIZE_PER_COLUMN_PARTITION_PER_SAMPLE = 16*1024L;
  private final long DEFAULT_SEGMENT_SIZE = 1048576L;
  private final String DEFAULT_ARRAY_NAME_PREFIX = "genomicsdb_array";

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_WORKSPACE_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_WORKSPACE_SHORT_NAME,
    doc = "Workspace where the database will be persisted. " +
      "Must be either a local path or full path of a " +
      "parallel file system like Lustre or NFS")
  private String workspace = "";

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_ARRAY_LONG_NAME,
  shortName = StandardArgumentDefinitions.GENOMICSDB_ARRAY_SHORT_NAME,
  doc = "TileDB array name used by GenomicsDB",
  optional = true)
  private String arrayName = DEFAULT_ARRAY_NAME_PREFIX;

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_SEGMENT_SIZE_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_SEGMENT_SIZE_SHORT_NAME,
    doc = "Buffer size in bytes allocated for GenomicsDB attributes during " +
      "import. Should be large enough to hold one variant data for one sample",
    optional = true)
  private
  long segmentSize = DEFAULT_SEGMENT_SIZE;

  @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
  shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
  doc = "GVCF files to be imported to GenomicsDB")
  private List<String> variantFileNames;

  @Argument(fullName = StandardArgumentDefinitions.SIZE_PER_COLUMN_PARTITION_PER_SAMPLE_LONG_NAME,
  shortName = StandardArgumentDefinitions.SIZE_PER_COLUMN_PARTITION_PER_SAMPLE_SHORT_NAME,
  doc = "Buffer size in bytes allocated to GenomicsDB variant readers. Default is 16KBytes." +
    " Larger values are better as smaller values cause frequent disk writes",
  optional = true)
  private long sizePerColumnPartitionPerSample = DEFAULT_SIZE_PER_COLUMN_PARTITION_PER_SAMPLE;

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_VID_MAP_FILE_LONG_NAME,
  shortName = StandardArgumentDefinitions.GENOMICSDB_VID_MAP_FILE_SHORT_NAME,
  doc = "JSON file including INFO/FORMAT/FILTER fields from headers and contig maps")
  private String vidMapJSONFile = "";

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_CALLSET_MAP_FILE_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_CALLSET_MAP_FILE_SHORT_NAME,
    doc = "JSON file including callset to GenomicsDB row maps")
  private String callsetMapJSONFile = "";

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_CREATE_WORKSPACE_LONG_NAME,
  shortName = StandardArgumentDefinitions.GENOMICSDB_CREATE_WORKSPACE_SHORT_NAME,
  doc = "Boolean flag - true means it will create GenomicsDB workspace, otherwise not",
  optional = true)
  private Boolean createWorkspace = false;

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

  /**
   * From the command line validate
   * 1. GenomicsDB workspace exists
   * 2. GenomicsDB workspace is a valid one
   *
   * @return  Error messages for the exceptions
   */
  @Override
  protected String[] customCommandLineValidation() {
    String[] messages = new String[2];
    int index = 0;

    if (workspace.isEmpty()) {
      messages[index] = " Must specify a workspace with -GW or --GenomicsDBWorkspace option";
      return messages;
    }

    File workspaceDir = new File(workspace);

    if (!createWorkspace) {
      if (!Files.isDirectory(workspaceDir.toPath()))
        messages[index++] = "No such directory or workspace found: " + workspace;

      boolean check = new File(workspace + "/__tiledb_workspace.tdb").exists();

      if (Files.isDirectory(workspaceDir.toPath()) && !check) {
        messages[index++] = "Not a valid GenomicsDB workspace: " + workspace;
      }
    }
    
    if (sizePerColumnPartitionPerSample < 1024L) {
      messages[index++] = "Size per column partition per sample is tool small." +
        " Consider using larger size than 10KBytes";
    }

    return index==0 ? null : messages;
  }

  /**
   * Before traversal, fix configuration parameters and initialize
   * GenomicsDB. Hard-coded to handle only VCF files and headers
   */
  @Override
  public void onTraversalStart() {

    File workspaceDir = new File(workspace);

    if (createWorkspace) {

      int ret = GenomicsDBImporter.createTileDBWorkspace(workspaceDir.getAbsolutePath());

      if (ret > 0) {
        try {
          throw new IOException("Directory " + workspaceDir + " already exists");
        } catch (IOException e) {
          e.printStackTrace();
          return;
        }
      } else if (ret < 0) {
        try {
          throw new IOException("Error creating GenomicsDB workspace");
        } catch (IOException e) {
          e.printStackTrace();
          return;
        }
      }
    }

    initializeIntervals();

    List<VCFHeader> headers = new ArrayList<>();

    Map<String, FeatureReader<VariantContext>> sampleToVCFMap = new HashMap<>();
    for (String variantFile : variantFileNames) {
      String variantFileAbsPath = new File(variantFile).getAbsolutePath();
      AbstractFeatureReader<VariantContext, LineIterator> reader =
        AbstractFeatureReader.getFeatureReader(variantFileAbsPath, new VCFCodec(), true);
      headers.add((VCFHeader)reader.getHeader());

      // We assume only one sample per file
      String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
      sampleToVCFMap.put(sampleName, reader);
    }

    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);
    sizePerColumnPartitionPerSample *= sampleToVCFMap.size();

    logger.info("Writing data to array - " + workspace + "/" + arrayName);

    try {
      importer = new GenomicsDBImporter(sampleToVCFMap, mergedHeader,
        intervals.get(0), workspace, arrayName,
        sizePerColumnPartitionPerSample, segmentSize,
        vidMapJSONFile, callsetMapJSONFile);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
   * from scratch can extend this class directly and implement this method. Walker authors should
   * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
   * Walker base classes implement the various kinds of traversals for you.
   */
  @Override
  public void traverse() {
    try {
      importer.importBatch();
      logger.info("Current batch written to disk");
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Loads our intervals using the best available sequence
   * dictionary (as returned by {@link #getBestAvailableSequenceDictionary})
   * to parse/verify them. Does nothing if no intervals were specified.
   */
  private void initializeIntervals() {
    if ( intervalArgumentCollection.intervalsSpecified() ) {
      final SAMSequenceDictionary intervalDictionary = getBestAvailableSequenceDictionary();
      if ( intervalDictionary == null ) {
        throw new UserException("We require at least one input source that " +
          "has a sequence dictionary (reference or reads) when intervals are specified");
      }

      intervals = new ArrayList<>();

      List<SimpleInterval> simpleIntervalList =
        intervalArgumentCollection.getIntervals(intervalDictionary);

      for (SimpleInterval simpleInterval : simpleIntervalList) {
        intervals.add(new ChromosomeInterval(simpleInterval.getContig(),
          simpleInterval.getStart(), simpleInterval.getEnd()));
      }
    }
  }

  /**
   * Close the GenomicsDB importer
   * once all variants are written
   */
  @Override
  public void closeTool() {
    assert importer.isDone();
    logger.info("Successfully imported " + variantFileNames.size() + " callsets");
  }
}

