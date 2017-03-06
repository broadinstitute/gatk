package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.ChromosomeInterval;
import com.intel.genomicsdb.GenomicsDBImportConfiguration;
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
import org.apache.commons.collections.OrderedBidiMap;
import org.apache.commons.collections.map.HashedMap;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ColumnPartitionArgumentCollection;
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

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_WORKSPACE_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_WORKSPACE_SHORT_NAME,
    doc = "Workspace where the database will be persisted. " +
      "Must be either a local path or full path of a " +
      "parallel file system like Lustre or NFS",
    optional = false)
  private String workspace = "";

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_PRODUCE_COMBINED_VCF_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_PRODUCE_COMBINED_VCF_SHORT_NAME,
    doc = "Produce Combined VCF using GenomicsDB import Tool")
  Boolean produceCombinedVCF = false;

  @ArgumentCollection
  ColumnPartitionArgumentCollection columnPartitionArgumentCollection =
    new ColumnPartitionArgumentCollection();

  @Argument(fullName = StandardArgumentDefinitions.COLUMN_PARTITION_BUFFER_LIMIT_LONG_NAME,
    shortName = StandardArgumentDefinitions.COLUMN_PARTITION_BUFFER_LIMIT_SHORT_NAME,
    doc = "To produce a column major array, the program allocates abuffer for every" +
      " input sample/callSet.")
  int columnPartitionBuffer = 10485760;

  @Argument(fullName = StandardArgumentDefinitions.GENOMICSDB_SEGMENT_SIZE_LONG_NAME,
    shortName = StandardArgumentDefinitions.GENOMICSDB_SEGMENT_SIZE_SHORT_NAME,
    doc = "Buffer size in bytes allocated for TileDB attributes during the" +
      "loading process. Should be large enough to hold one cell worth of data.",
    optional = true)
  int segmentSize = 10485760;

  @Argument(fullName="variants", shortName="V",doc="variant files")
  List<FeatureInput<VariantContext>> variantContexts;

  @Override
  public boolean requiresIntervals() { return true; }

  @Override
    public boolean requiresReads() { return false; }

  @Override
  public boolean requiresReference() {
    return false;
  }

  GenomicsDBImporter importer;
  private boolean done = false;
  private List<SimpleInterval> intervalsForTraversal;

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
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

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

    File workspaceDir = new File(workspace);

    if (!Files.isDirectory(workspaceDir.toPath()))
      messages[index++] = "No such directory or workspace found: " + workspace;

    boolean check = new File(workspace + "/__tiledb_workspace.tdb").exists();
    if (Files.isDirectory(workspaceDir.toPath()) && !check) {
      messages[index++] = "Not a valid GenomicsDB workspace: " + workspace;
    }

    return index==0 ? null : messages;
  }

  /**
   * Before traversal, fix configuration parameters and initialize
   * GenomicsDB. Hard-coded to handle only VCF files and headers
   */
  @Override
  public void onTraversalStart() {

    List<VCFHeader> headers = new ArrayList<>();
    FeatureCodec<VariantContext,LineIterator> codec = new VCFCodec();

    Map<String, FeatureReader<VariantContext>> sampleToVCMap =
      new HashMap<String, FeatureReader<VariantContext>>();
    for (FeatureInput<VariantContext> variant : variantContexts) {
      AbstractFeatureReader<VariantContext, LineIterator> reader =
        AbstractFeatureReader.getFeatureReader(String.valueOf(variant), codec, false);
      headers.add((VCFHeader)reader.getHeader());

      // We assume only one sample per file
      String sampleName = ((VCFHeader) reader.getHeader()).getGenotypeSamples().get(0);
      sampleToVCMap.put(sampleName, reader);
    }

    Set<VCFHeaderLine> mergedHeader = VCFUtils.smartMergeHeaders(headers, true);

    List<GenomicsDBImportConfiguration.Partition> partitions = new ArrayList<>();
    GenomicsDBImportConfiguration.Partition.Builder pB =
      GenomicsDBImportConfiguration.Partition.newBuilder();
    GenomicsDBImportConfiguration.Partition p0 =
      pB
        .setBegin(0)
        .setTiledbWorkspace(workspace)
        .setTiledbArrayName("array0")
        .build();
    partitions.add(p0);

    GenomicsDBImportConfiguration.ImportConfiguration.Builder importBuilder =
      GenomicsDBImportConfiguration.ImportConfiguration.newBuilder();
    GenomicsDBImportConfiguration.ImportConfiguration importConfiguration =
      importBuilder
        .addAllColumnPartitions(partitions)
        .setCompressTiledbArray(true)
        .setSizePerColumnPartition(1000)
        .build();

    try {
      importer = new GenomicsDBImporter(
        sampleToVCMap,
        mergedHeader,
        toChromosomeInterval(intervalsForTraversal.get(0)),
        importConfiguration);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private ChromosomeInterval toChromosomeInterval(SimpleInterval interval) {
    return new ChromosomeInterval(
      interval.getContig(), interval.getStart(), interval.getEnd());
  }

  /**
   * Close the GenomicsDB importer
   * once all variants are written
   */
  @Override
  public void closeTool() {
    done = importer.isDone();
  }

  public void initializeIntervals() {
    if (intervalArgumentCollection.intervalsSpecified()) {
      final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
      if ( sequenceDictionary == null ) {
        throw new UserException("Sequence Dictionary is null");
      }

      intervalsForTraversal = intervalArgumentCollection.getIntervals(sequenceDictionary);
    }
  }
}

