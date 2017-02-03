package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.barclay.argparser.Argument;

public abstract class GenomicsDBCommandLineProgram extends CommandLineProgram {

  public static final String LOADER_JSON_FULL_NAME = "loaderJSONFile";
  public static final String LOADER_JSON_SHORT_NAME = "l";
  public static final String STREAM_ID_JSON_FULL_NAME = "streamIdJSONFile";
  public static final String STREAM_ID_JSON_SHORT_NAME = "s";
  public static final String BUFFER_CAPACITY_FULL_NAME = "bufferCapacity";
  public static final String BUFFER_CAPACITY_SHORT_NAME = "b";
  public static final String PARTITION_INDEX_FULL_NAME = "partitionIndex";
  public static final String PARTITION_INDEX_SHORT_NAME = "p";

  @Argument(fullName = BUFFER_CAPACITY_FULL_NAME,
    shortName = BUFFER_CAPACITY_SHORT_NAME,
    doc = "Buffer capacity for GenomicsDB update",
    optional = true)
  public Long bufferCapacity = 1024L;

  @Override
  public Object instanceMain(final String[] argv) {
    // First, we parse the commandline arguments, then we set important statics like VALIDATION_STRINGENCY, and
    // finally, we call into the normal instance main (post arg-parsing). If don't start with the argument parsing
    // we always get default values for VALIDATION_STRINGENCY, COMPRESSION_LEVEL, etc.
    if (!parseArgs(argv)) {
      //an information only argument like help or version was specified, just exit
      return 0;
    }

    // defer to parent to finish the initialization and starting the program.
    return instanceMainPostParseArgs();
  }
}
