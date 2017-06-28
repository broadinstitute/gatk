package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.fs.ContentSummary;
import org.apache.hadoop.fs.FileSystem;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary        = "Program to call SGA to perform local assembly and return assembled contigs if successful, " +
                         "or runtime error messages if the process erred for some breakpoints.",
        oneLineSummary = "Perform SGA-based local assembly on fasta files on Spark.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class RunSGAViaProcessBuilderOnSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc       = "Absolute path to SGA installation.",
              shortName = "sgaPath",
              fullName  = "fullPathToSGA",
              optional  = false)
    public String pathToSGA = null;

    @Argument(doc       = "An URI to the directory where all interleaved FASTQ files for putative breakpoints are located.",
              shortName = "in",
              fullName  = "inDir",
              optional  = false)
    public String pathToAllInterleavedFASTQFiles = null;

    @Argument(doc       = "A substring in the FASTQ file names that needs to be stripped out for retrieving the breakpoint ID",
              shortName = "scrub",
              fullName  = "subStringToStrip",
              optional  = false)
    public String subStringToStrip = null;

    @Argument(doc       = "An URI (prefix) to a directory to write results to. " +
                          "Breakpoints where local assembly were successful are saved in a directory prefix_0," +
                          "breakpoints where local assembly failed are saved in directory prefix_1",
              shortName = "out",
              fullName  = "outDirPrefix",
              optional  = false)
    public String outDirPrefix = null;

    @Argument(doc       = "To run k-mer based read correction with default parameters (as provided by SGA) or not.",
              shortName = "correct",
              fullName  = "kmerCorrectReads",
              optional  = true)
    public boolean runCorrection = false;

    @Argument(doc       = "If true, capture the stdout and stderr of the SGA processes along with the assembly output (valuable for debugging but hurts performance)",
              shortName = "captureStdIO",
              fullName  = "captureSGAStdOutAndError",
              optional  = true)
    public boolean enableSTDIOCapture = false;

    // a few hard-coded parameters for use in various SGA modules based on some tuning experiences.
    // subject to future changes as we see more test cases
    @VisibleForTesting static final int MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE = 55;
    @VisibleForTesting static final int FILTER_STEP_KMER_FREQUENCY_THREASHOLD = 1;
    @VisibleForTesting static final double OVERLAP_STEP_ERROR_RATE = 0.05;
    @VisibleForTesting static final double ASSEMBLE_STEP_MAX_GAP_DIVERGENCE = 0.0;
    @VisibleForTesting static final int ASSEMBLE_STEP_MIN_BRANCH_TAIL_LENGTH = 50;
    @VisibleForTesting static final int CUT_OFF_TERMINAL_BRANCHES_IN_N_ROUNDS = 0;
    @VisibleForTesting static final boolean TURN_ON_TRANSITIVE_REDUCTION = false;


    // for developer performance debugging use
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
    private static final Logger logger = LogManager.getLogger(RunSGAViaProcessBuilderOnSpark.class);

    @Override
    public void runTool(final JavaSparkContext ctx){

        logger.debug("SGAOnSpark_debug: Start job at " + dateFormat.format(new Date()));

        // first load RDD of pair that has path to FASTQ file path as its first and FASTQ file contents as its second
        final JavaPairRDD<String, String> fastqContentsForEachBreakpoint = loadFASTQFiles(ctx, pathToAllInterleavedFASTQFiles);

        final JavaPairRDD<Long, SGAAssemblyResult> assembly = fastqContentsForEachBreakpoint.mapToPair(entry -> performAssembly(entry, subStringToStrip, pathToSGA, runCorrection, enableSTDIOCapture));

        validateAndSaveResults(assembly, outDirPrefix);

        logger.debug("SGAOnSpark_debug: Finish job at " + dateFormat.format(new Date()));
    }

    /**
     * Load the FASTQ files in the user specified directory and returns an RDD that satisfies the same requirement
     * as described in {@link JavaSparkContext#wholeTextFiles(String, int)}.
     * @param ctx
     * @param pathToAllInterleavedFASTQFiles path to the directory where all FASTQ files to perform local assembly upon are located
     * @throws GATKException when getting the file count in the specified directory
     * @return
     */
    private static JavaPairRDD<String, String> loadFASTQFiles(final JavaSparkContext ctx, final String pathToAllInterleavedFASTQFiles){
        try{
            final FileSystem hadoopFileSystem = FileSystem.get(ctx.hadoopConfiguration());
            final ContentSummary cs = hadoopFileSystem.getContentSummary(new org.apache.hadoop.fs.Path(pathToAllInterleavedFASTQFiles));
            final int fileCount = (int) cs.getFileCount();
            return ctx.wholeTextFiles(pathToAllInterleavedFASTQFiles, fileCount);
        }catch (final IOException e){
            throw new GATKException(e.getMessage());
        }
    }

    /**
     * Validates the returned result from running the local assembly pipeline:
     *   if all steps executed successfully, the contig file is nonnull so we save the contig file and discard the runtime information
     *   if any SGA step returns non-zero code, the contig file is null so we save the runtime information for that break point
     * @param results       the local assembly result and its associated breakpoint ID
     * @param outputDir     output directory to save the contigs (if assembly succeeded) or runtime info (if erred)
     */
    private static void validateAndSaveResults(final JavaPairRDD<Long, SGAAssemblyResult> results, final String outputDir){

        final JavaPairRDD<Long, SGAAssemblyResult> cachedResults = results.cache(); // cache because Spark doesn't have an efficient RDD.split(predicate) yet

        // save fasta file contents or failure message
        final JavaPairRDD<Long, SGAAssemblyResult> success = cachedResults.filter(entry -> entry._2().assembledContigs!=null);
        final JavaPairRDD<Long, SGAAssemblyResult> failure = cachedResults.filter(entry -> entry._2().assembledContigs==null);

        if(!success.isEmpty()){
            success.map(entry -> entry._1().toString() + "\t" + entry._2().assembledContigs.toPackedFasta())
                    .saveAsTextFile(outputDir+"_0");
        }

        if(!failure.isEmpty()){
            failure.map(entry ->  entry._1().toString() + "\n" + entry._2().collectiveRuntimeInfo.toString())
                    .saveAsTextFile(outputDir+"_1");
        }
    }

    /**
     * Performs assembly on the FASTA files pointed to by the URI that is associated with the breakpoint identified by the long ID.
     * Actual assembly work is delegated to other functions.
     * @param fastqOfABreakpoint    the (partial) URI to the FASTQ file and FASTQ file contents as String
     * @param subStringInFilenameToScrub the part in a file name that must be stripped out to extract breakpoint ID, e.g. "assembly1234" -> 1234
     * @param sgaPath               full path to SGA
     * @param runCorrections        user's decision to run SGA's corrections (with default parameter values) or not
     * @param enableSTDIOCapture    to enable capture of stderr & stdout of SGA processes or not
     * @return                      failure message (if process erred) or contig FASTA file contents (if process succeeded) associated with the breakpoint ID
     * @throws IOException          if fails to create temporary directory on local filesystem or fails to parse contig FASTA file
     */
    @VisibleForTesting
    static Tuple2<Long, SGAAssemblyResult> performAssembly(final Tuple2<String, String> fastqOfABreakpoint,
                                                           final String subStringInFilenameToScrub,
                                                           final String sgaPath,
                                                           final boolean runCorrections,
                                                           final boolean enableSTDIOCapture)
            throws IOException{

        final Tuple2<Long, File> localFASTQFileForOneBreakpoint = writeToLocal(fastqOfABreakpoint, subStringInFilenameToScrub);

        final SGAAssemblyResult assembledContigsFileAndRuntimeInfo = runSGAModulesInSerial(sgaPath, localFASTQFileForOneBreakpoint._2(), runCorrections, enableSTDIOCapture);

        return new Tuple2<>(localFASTQFileForOneBreakpoint._1(), assembledContigsFileAndRuntimeInfo);
    }

    /**
     * Utility function that unloads the FASTQ contents for a breakpoint to a local file for later consumption by SGA.
     * @param oneBreakPoint input for one breakpoint, where the first is the path to the FASTQ file and the second is the FASTQ file's content
     * @return              the breakpoint ID and with the FASTQ file contents dumped to a local File
     * @throws IOException  if fails to create the temporary directory or fails to write to local file
     */
    @VisibleForTesting
    static Tuple2<Long, File> writeToLocal(final Tuple2<String, String> oneBreakPoint, final String subStringToStripout) throws IOException{

        final String fastqFilename = FilenameUtils.getName(oneBreakPoint._1());

        final File localTempWorkingDir = Files.createTempDirectory( fastqFilename + "_" ).toAbsolutePath().toFile();
        localTempWorkingDir.deleteOnExit();

        final File localFASTQFile =  new File(localTempWorkingDir, fastqFilename);
        FileUtils.writeStringToFile(localFASTQFile, oneBreakPoint._2());

        final Long breakpointID = Long.parseLong(FilenameUtils.getBaseName(oneBreakPoint._1()).replace(subStringToStripout, ""));

        return new Tuple2<>(breakpointID, localFASTQFile);
    }

    /**
     * Linear pipeline for running the SGA local assembly process on a particular FASTQ file for its associated putative breakpoint.
     *
     * @param sgaPathString     full path to the SGA program
     * @param rawFASTQFile      local FASTQ file living in a temp local working dir
     * @param runCorrections    to run SGA correction steps--correct, filter, rmdup, merge--or not
     * @param enableSTDIOCapture to enable capture of stderr & stdout of SGA processes or not
     * @return                  the result accumulated through running the pipeline, where the contigs file name could be null if the process erred.
     * @throws IOException      if fails to parse the final contigs file
     */
    @VisibleForTesting
    static SGAAssemblyResult runSGAModulesInSerial(final String sgaPathString,
                                                   final File rawFASTQFile,
                                                   final boolean runCorrections,
                                                   final boolean enableSTDIOCapture)
            throws IOException{

        logger.debug("SGAOnSpark_debug: start actual assembly process for " + rawFASTQFile.getName() + " at " + dateFormat.format(new Date()));

        final Path sgaPath = Paths.get(sgaPathString);

        final File workingDir = rawFASTQFile.getParentFile();

        // the index module is used frequently, so make single instance and pass around
        final SGAModule indexer = new SGAModule("index");
        final List<String> indexerArgs = new ArrayList<>();
        indexerArgs.add("--algorithm"); indexerArgs.add("ropebwt");
        indexerArgs.add("--check");
        indexerArgs.add("");

        // collect runtime information along the way, if user requests
        final List<SGAModule.RuntimeInfo> runtimeInfo = new ArrayList<>();

        final Map<String, String> moduleArgsAndValues = new LinkedHashMap<>(10); // none of the SGA modules used in the following have more than 10 sensibly tunable parameters that affects assembly results

        String nameOfFileToFeedToNextStep = null;
        {//preprocess
            moduleArgsAndValues.clear();
            moduleArgsAndValues.put("--pe-mode", "2");
            moduleArgsAndValues.put("--out", FilenameUtils.getBaseName(rawFASTQFile.getName()) + ".pp.fa");

            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "preprocess", ".pp.fa", rawFASTQFile, workingDir, indexer, indexerArgs, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        if(runCorrections){//optional kmer based error correction
            moduleArgsAndValues.clear();
            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "correct", ".ec.fa", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, indexer, indexerArgs, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        {//filter
            moduleArgsAndValues.clear();
            moduleArgsAndValues.put("--kmer-threshold", String.valueOf(FILTER_STEP_KMER_FREQUENCY_THREASHOLD));
            moduleArgsAndValues.put("--homopolymer-check", "");
            moduleArgsAndValues.put("--low-complexity-check", "");

            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "filter", ".filter.pass.fa", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, null, null, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        {//merge
            moduleArgsAndValues.clear();
            moduleArgsAndValues.put("--min-overlap", String.valueOf(MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));

            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "fm-merge", ".merged.fa", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, indexer, indexerArgs, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        {//deduplicate
            moduleArgsAndValues.clear();
            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "rmdup", ".rmdup.fa", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, indexer, indexerArgs, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        {//overlap
            moduleArgsAndValues.clear();
            moduleArgsAndValues.put("--min-overlap", String.valueOf(MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));
            moduleArgsAndValues.put("--error-rate", String.valueOf(OVERLAP_STEP_ERROR_RATE));

            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "overlap", ".asqg.gz", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, null, null, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){ return new SGAAssemblyResult(null, runtimeInfo); }
        }

        {//assembly
            moduleArgsAndValues.clear();
            moduleArgsAndValues.put("--cut-terminal",          String.valueOf(CUT_OFF_TERMINAL_BRANCHES_IN_N_ROUNDS));
            moduleArgsAndValues.put("--max-gap-divergence",    String.valueOf(ASSEMBLE_STEP_MAX_GAP_DIVERGENCE));
            moduleArgsAndValues.put("--min-overlap",           String.valueOf(MIN_OVERLAP_IN_FILTER_OVERLAP_ASSEMBLE));
            moduleArgsAndValues.put("--min-branch-length",     String.valueOf(ASSEMBLE_STEP_MIN_BRANCH_TAIL_LENGTH));
            moduleArgsAndValues.put("--out-prefix",            nameOfFileToFeedToNextStep.replace(".asqg.gz", ""));

            nameOfFileToFeedToNextStep = runSGAModule(sgaPath, "assemble", "-contigs.fa", new File(workingDir, nameOfFileToFeedToNextStep), workingDir, null, null, turnCmdLineArgsKeyValuePairIntoList(moduleArgsAndValues), runtimeInfo, enableSTDIOCapture);
            if( null == nameOfFileToFeedToNextStep ){
                return new SGAAssemblyResult(null, runtimeInfo);
            }else{
                nameOfFileToFeedToNextStep = nameOfFileToFeedToNextStep.replace(".asqg", ""); // strip out substring ".asqg"
            }
        }

        logger.debug("SGAOnSpark_debug: finished actual assembly process for " + rawFASTQFile.getName() + " at " + dateFormat.format(new Date()));

        // if code reaches here, all steps in the SGA pipeline went smoothly,
        // but the following conversion from File to ContigsCollection may still err
        final File assembledContigsFile = new File(workingDir, nameOfFileToFeedToNextStep);
        final List<String> contigsFASTAContents = !assembledContigsFile.exists() ? null : Files.readAllLines(Paths.get(assembledContigsFile.getAbsolutePath()));
        if(null==contigsFASTAContents){
            return new SGAAssemblyResult(null, runtimeInfo);
        }else{
            return new SGAAssemblyResult(contigsFASTAContents, runtimeInfo);
        }
    }

    /**
     * Call the right SGA module, log runtime information, and return the output file name if succeed.
     * If process erred, the string returned is null.
     * @param moduleName            SGA module name to be run
     * @param inputFile             file to be fed to SGA module
     * @param sgaPath               full path to the SGA program
     * @param workingDir            directory the SGA pipeline is working in
     * @param indexer               module representing SGA index
     * @param indexerArgs           arguments used by SGA index
     * @param collectedRuntimeInfo  runtime information collected along the process
     * @param enableSTDIOCapture    to enable capture of stdout & stderr of SGA processes or not
     * @return                      the name of the file produced by running this SGA module, {@code null} if module failed to run
     */
    @VisibleForTesting
    static String runSGAModule(final Path sgaPath,
                               final String moduleName,
                               final String extensionToAppend,
                               final File inputFile,
                               final File workingDir,
                               final SGAModule indexer,
                               final List<String> indexerArgs,
                               final List<String> moduleArgs,
                               final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                               final boolean enableSTDIOCapture){

        final SGAModule module = new SGAModule(moduleName);

        moduleArgs.add(inputFile.getName());

        final SGAModule.RuntimeInfo moduleInfo = module.run(sgaPath, workingDir, moduleArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(moduleInfo);

        final String outputFileName = FilenameUtils.getBaseName(inputFile.getName()) + extensionToAppend;

        if(null!=indexer){
            indexerArgs.set(indexerArgs.size()-1, outputFileName);
            final SGAModule.RuntimeInfo indexerInfo = indexer.run(sgaPath, workingDir, indexerArgs, enableSTDIOCapture);
            collectedRuntimeInfo.add(indexerInfo);
        }

        // check module return status
        final SGAModule.RuntimeInfo.ReturnStatus returnStatus = collectedRuntimeInfo.get(collectedRuntimeInfo.size()-1).returnStatus;
        if(!(returnStatus.equals(SGAModule.RuntimeInfo.ReturnStatus.SUCCESS))){
            return null;
        }else{
            return outputFileName;
        }
    }

    @VisibleForTesting
    static List<String> turnCmdLineArgsKeyValuePairIntoList(final Map<String, String> keyValuePairs){

        return keyValuePairs.entrySet().stream().flatMap(p -> {final String val = p.getValue();
                                                               return val.isEmpty() ? Stream.of(p.getKey()) : Stream.of(p.getKey(), val);
                                                               })
                                                .collect(Collectors.toList());
    }

    /**
     * Final return type of the whole process of SGA local assembly.
     * assembledContigFile is the file containing the assembled contigs, if the process executed successfully, or null if not.
     * runtimeInformation contains the runtime information logged along the process up until the process erred, if errors happen,
     *   or until the last step if no errors occur along the line.
     *   The list is organized along the process of executing the assembly pipeline.
     */
    @VisibleForTesting
    static final class SGAAssemblyResult implements Serializable{
        private static final long serialVersionUID = 1L;

        public final ContigsCollection assembledContigs;
        public final List<SGAModule.RuntimeInfo> collectiveRuntimeInfo;

        public SGAAssemblyResult(final List<String> fastaContents, final List<SGAModule.RuntimeInfo> collectedRuntimeInfo){
            this.assembledContigs      =  null== fastaContents ? null : new ContigsCollection(fastaContents);
            this.collectiveRuntimeInfo = collectedRuntimeInfo;
        }

        @Override
        public String toString(){
            if(null!=assembledContigs){
                return assembledContigs.toString();
            }else{
                return StringUtils.join(collectiveRuntimeInfo.stream().map(info -> info.toString()), "\n");
            }
        }
    }

}