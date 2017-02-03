package org.broadinstitute.hellbender.tools.genomicsdb;

import com.intel.genomicsdb.VCF2TileDB;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.GenomicsDBCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.programgroups.GenomicsDBProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.json.simple.parser.ContainerFactory;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.IOException;
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
  summary = "Import VCF to GenomicsDB",
  oneLineSummary = "Import VCF to GenomicsDB",
  programGroup = GenomicsDBProgramGroup.class
)
public final class GenomicsDBImport extends GenomicsDBCommandLineProgram {

  @Argument(doc = "Input Loader JSON File. This configuration file contains execution parameters " +
    "for loading variants to GenomicsDB. Please refer to " +
    "https://github.com/Intel-HLS/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB#execution-parameters-for-the-import-program" +
    "for more information",
    fullName = LOADER_JSON_FULL_NAME,
    shortName = LOADER_JSON_SHORT_NAME,
    common=true)
  private String loaderJSONFile = "";

  @Argument(doc = "Input JSON File with Stream Ids. This configuration file contains stream names " +
    "and corresponding file names for all samples",
    fullName = STREAM_ID_JSON_FULL_NAME,
    shortName = STREAM_ID_JSON_SHORT_NAME,
    common=true)
  private String streamIdJSONFile = "";

  @Argument(doc = "Rank of the process",
    fullName = PARTITION_INDEX_FULL_NAME,
    shortName = PARTITION_INDEX_SHORT_NAME,
    optional = true)
  private int partitionIndex = 0;

  /**
   * Do the work after command line has been parsed. RuntimeException may be
   * thrown by this method, and are reported appropriately.
   *
   * @return the return value or null is there is none.
   */
  @Override
  protected Object doWork() {

    VCF2TileDB loader = new VCF2TileDB(loaderJSONFile, partitionIndex);
    try {
      JSONParser parser = new JSONParser();
      Object obj = parser.parse(new FileReader(streamIdJSONFile), new LinkedHashFactory<String,String>());
      LinkedHashMap<?, ?> streamNameToFileName = (LinkedHashMap<?,?>) obj;
      long rowIdx = 0;
      for (Object currObj : streamNameToFileName.entrySet()) {
        Map.Entry<?, ?> entry = (Map.Entry<?, ?>) currObj;
        VCFFileStreamInfo currInfo = new VCFFileStreamInfo(entry.getValue().toString(),
          loaderJSONFile, partitionIndex);
        LinkedHashMap<Integer, com.intel.genomicsdb.VCF2TileDB.SampleInfo> sampleIndexToInfo =
          new LinkedHashMap<>();

        // Increments row index after every 
        rowIdx = com.intel.genomicsdb.VCF2TileDB.initializeSampleInfoMapFromHeader(sampleIndexToInfo,
          currInfo.mVCFHeader, rowIdx);
        loader.addSortedVariantContextIterator(entry.getKey().toString(),
          currInfo.mVCFHeader, currInfo.mIterator,
          bufferCapacity, VariantContextWriterBuilder.OutputType.BCF_STREAM,
          sampleIndexToInfo); //pass sorted VC iterators
      }
      loader.importBatch();
      assert loader.isDone();
    } catch (IOException e) {
      throw new UserException(streamIdJSONFile + ": No such file or directory" + e);
    } catch (ParseException pe)  {
      throw new UserException("Parse error for file: " + streamIdJSONFile + pe);
    }
    return null;
  }

  @Override
  protected String[] customCommandLineValidation() {
    String[] messages = new String[2];
    int index = 0;
    
    if (loaderJSONFile.isEmpty()) {
      messages[index++] = "Must provide an input loader JSON file. " +
        "Please visit https://github.com/Intel-HLS/GenomicsDB/wiki/Java-interface-for-importing-VCF-CSV-files-into-TileDB-GenomicsDB" +
        " for more information";
    }

    if (streamIdJSONFile.isEmpty()) {
      messages[index++] = "Must provide a stream JSON file. ";
    }
    
    if (index == 0) return null;
    return messages;
  }
}


/**
 * Factory object to maintain order of keys in simple JSON parsing
 * - use LinkedHashMap
 */
class LinkedHashFactory<E, S> implements ContainerFactory
{
  @Override
  public List<E> creatArrayContainer()
  {
    return new ArrayList<>();
  }

  @Override
  public Map<E,S> createObjectContainer()
  {
    return new LinkedHashMap<>();
  }
}

/**
 * A container class to maintain the input VCF streams
 */
class VCFFileStreamInfo
{
  VCFHeader mVCFHeader = null;
  Iterator<VariantContext> mIterator = null;

  /**
   * Constructor
   * @param fileName path to VCF file
   */
  VCFFileStreamInfo(final String fileName,
                    final String loaderJSONFile,
                    final int partitionIndex) throws IOException, ParseException
  {
    AbstractFeatureReader<VariantContext, LineIterator> reader =
      AbstractFeatureReader.getFeatureReader(fileName, new VCFCodec(), false);
    mVCFHeader = (VCFHeader)(reader.getHeader());
    mIterator = com.intel.genomicsdb.VCF2TileDB.columnPartitionIterator(
      reader, loaderJSONFile, partitionIndex);
  }
}