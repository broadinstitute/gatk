package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.hadoop.fs.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import scala.Tuple2;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

@CommandLineProgramProperties(
        summary        = "Program to call SGA to perform local assembly and return assembled contigs if successful, " +
                         "or runtime error messages if the process erred for some breakpoints.",
        oneLineSummary = "Perform SGA-based local assembly on fasta files on Spark.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
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
     *   if any sga step returns non-zero code, the contig file is null so we save the runtime information for that break point
     * @param results       the local assembly result and its associated breakpoint ID
     * @param outputDir     output directory to save the contigs (if assembly succeeded) or runtime info (if erred)
     */
    private static void validateAndSaveResults(final JavaPairRDD<Long, SGAAssemblyResult> results, final String outputDir){

        final JavaPairRDD<Long, SGAAssemblyResult> cachedResults = results.cache(); // cache because Spark doesn't have an efficient RDD.split(predicate) yet

        // save fasta file contents or failure message
        final JavaPairRDD<Long, SGAAssemblyResult> success = cachedResults.filter(entry -> entry._2().assembledContigs!=null);
        final JavaPairRDD<Long, SGAAssemblyResult> failure = cachedResults.filter(entry -> entry._2().assembledContigs==null);

        if(!success.isEmpty()){
            success.map(entry -> entry._1().toString() + "\t" + entry._2().assembledContigs.toOneLineString())
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
     * @param sgaPath           full path to the SGA program
     * @param rawFASTQFile      local FASTQ file living in a temp local working dir
     * @param runCorrections    to run SGA correction steps--correct, filter, rmdup, merge--or not
     * @param enableSTDIOCapture to enable capture of stderr & stdout of SGA processes or not
     * @return                  the result accumulated through running the pipeline, where the contigs file name could be null if the process erred.
     * @throws IOException      if fails to parse the final contigs file
     */
    @VisibleForTesting
    static SGAAssemblyResult runSGAModulesInSerial(final String sgaPath,
                                                   final File rawFASTQFile,
                                                   final boolean runCorrections,
                                                   final boolean enableSTDIOCapture)
            throws IOException{

        logger.debug("SGAOnSpark_debug: start actual assembly process for " + rawFASTQFile.getName() + " at " + dateFormat.format(new Date()));

        final Path sga = Paths.get(sgaPath);

        final File tempWorkingDir = rawFASTQFile.getParentFile();

        // the index module is used frequently, so make single instance and pass around
        final SGAModule indexer = new SGAModule("index");
        final List<String> indexerArgs = new ArrayList<>();
        indexerArgs.add("--algorithm"); indexerArgs.add("ropebwt");
        indexerArgs.add("--check");
        indexerArgs.add("");

        // collect runtime information along the way, if user requests
        final List<SGAModule.RuntimeInfo> runtimeInfo = new ArrayList<>();

        String preppedFileName = runAndStopEarly("preprocess", rawFASTQFile, sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
        if( null == preppedFileName ){
            return new SGAAssemblyResult(null, runtimeInfo);
        }

        if(runCorrections){
            preppedFileName = runAndStopEarly("correct", new File(tempWorkingDir, preppedFileName), sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
            if( null == preppedFileName ){
                return new SGAAssemblyResult(null, runtimeInfo);
            }
        }

        final File preppedFile = new File(tempWorkingDir, preppedFileName);
        final String filteredFileName = runAndStopEarly("filter", preppedFile, sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
        if(null == filteredFileName){
            return new SGAAssemblyResult(null, runtimeInfo);
        }

        final File fileToMerge = new File(tempWorkingDir, filteredFileName);
        final String fileNameToDedup = runAndStopEarly("fm-merge", fileToMerge, sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
        if(null == fileNameToDedup){
            return new SGAAssemblyResult(null, runtimeInfo);
        }

        final File fileToDeduplicate = new File(tempWorkingDir, fileNameToDedup);
        final String fileNameToAssemble = runAndStopEarly("rmdup", fileToDeduplicate, sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
        if(null == fileNameToAssemble){
            return new SGAAssemblyResult(null, runtimeInfo);
        }

        final File fileToAssemble   = new File(tempWorkingDir, fileNameToAssemble);
        final String contigsFileName = runAndStopEarly("assemble", fileToAssemble, sga, tempWorkingDir, indexer, indexerArgs, runtimeInfo, enableSTDIOCapture);
        if(null == contigsFileName){
            return new SGAAssemblyResult(null, runtimeInfo);
        }

        logger.debug("SGAOnSpark_debug: finished actual assembly process for " + rawFASTQFile.getName() + " at " + dateFormat.format(new Date()));

        // if code reaches here, all steps in the SGA pipeline went smoothly,
        // but the following conversion from File to ContigsCollection may still err
        final File assembledContigsFile = new File(tempWorkingDir, contigsFileName);
        final List<String> contigsFASTAContents = (null==assembledContigsFile) ? null : Files.readAllLines(Paths.get(assembledContigsFile.getAbsolutePath()));
        if(null==contigsFASTAContents){
            return new SGAAssemblyResult(null, runtimeInfo);
        }else{
            return new SGAAssemblyResult(contigsFASTAContents, runtimeInfo);
        }
    }

    /**
     * Call the right sga module, log runtime information, and return the output file name if succeed.
     * If process erred, the string returned is null.
     * @param moduleName            SGA module name to be run
     * @param inputFASTQFile        FASTQ file tobe fed to SGA module
     * @param sgaPath               full path to the SGA program
     * @param workingDir            directory the SGA pipeline is working in
     * @param indexer               module representing SGA index
     * @param indexerArgs           arguments used by SGA index
     * @param collectedRuntimeInfo  runtime information collected along the process
     * @param enableSTDIOCapture    to enable capture of stdout & stderr of sag processes or not
     * @return                      the name of file produced by running this SGA module
     */
    private static String runAndStopEarly(final String moduleName,
                                          final File inputFASTQFile,
                                          final Path sgaPath,
                                          final File workingDir,
                                          final SGAModule indexer,
                                          final List<String> indexerArgs,
                                          final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                          final boolean enableSTDIOCapture){

        String filenameToReturn = null;
        if(moduleName.equalsIgnoreCase("preprocess")){
            filenameToReturn = runSGAPreprocess(sgaPath, inputFASTQFile, workingDir, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
        }else if(moduleName.equalsIgnoreCase("correct")){
            filenameToReturn = runSGACorrect(sgaPath, inputFASTQFile, workingDir, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
        }else if(moduleName.equalsIgnoreCase("filter")){
            filenameToReturn = runSGAFilter(sgaPath, inputFASTQFile, workingDir, collectedRuntimeInfo, enableSTDIOCapture);
        }else if(moduleName.equalsIgnoreCase("rmdup")){
            filenameToReturn = runSGARmDuplicate(sgaPath, inputFASTQFile, workingDir, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
        }else if(moduleName.equalsIgnoreCase("fm-merge")){
            filenameToReturn = runSGAFMMerge(sgaPath, inputFASTQFile, workingDir, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
        }else if(moduleName.equalsIgnoreCase("assemble")){
            filenameToReturn = runSGAOverlapAndAssemble(sgaPath, inputFASTQFile, workingDir, collectedRuntimeInfo, enableSTDIOCapture);
        }else{
            throw new GATKException("Wrong module called"); // should never occur, implementation mistake
        }

        final SGAModule.RuntimeInfo.ReturnStatus returnStatus = collectedRuntimeInfo.get(collectedRuntimeInfo.size()-1).returnStatus;

        if(!(returnStatus.equals( SGAModule.RuntimeInfo.ReturnStatus.SUCCESS))){
            return null;
        }else{
            return filenameToReturn;
        }
    }

    @VisibleForTesting
    static String runSGAPreprocess(final Path sgaPath,
                                   final File inputFASTQFile,
                                   final File outputDirectory,
                                   final SGAModule indexer,
                                   final List<String> indexerArgs,
                                   final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                   final boolean enableSTDIOCapture){

        final String prefix = FilenameUtils.getBaseName(inputFASTQFile.getName());
        final String preprocessedFASTAFileName = prefix+".pp.fa";

        final SGAModule preprocess = new SGAModule("preprocess");
        final List<String> ppArgs = new ArrayList<>();
        ppArgs.add("--pe-mode");    ppArgs.add("2");
        ppArgs.add("--pe-orphans"); ppArgs.add(prefix+".pp.orphan.fa");
        ppArgs.add("--out");        ppArgs.add(preprocessedFASTAFileName);
        ppArgs.add(inputFASTQFile.getName());

        final SGAModule.RuntimeInfo ppInfo = preprocess.run(sgaPath, outputDirectory, ppArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(ppInfo);

        indexerArgs.set(indexerArgs.size()-1, preprocessedFASTAFileName);
        final SGAModule.RuntimeInfo indexerInfo = indexer.run(sgaPath, outputDirectory, indexerArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(indexerInfo);

        return preprocessedFASTAFileName;
    }

    @VisibleForTesting
    static String runSGACorrect(final Path sgaPath,
                                final File inputFASTAFile,
                                final File outputDirectory,
                                final SGAModule indexer,
                                final List<String> indexerArgs,
                                final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                final boolean enableSTDIOCapture){
        return runSimpleModuleFollowedByIndexing(sgaPath, "correct", ".ec.fa", inputFASTAFile, outputDirectory, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
    }

    @VisibleForTesting
    static String runSGAFilter(final Path sgaPath,
                               final File inputFASTAFile,
                               final File outputDirectory,
                               final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                               final boolean enableSTDIOCapture){

        final String prefix = FilenameUtils.getBaseName(inputFASTAFile.getName());
        final SGAModule filter = new SGAModule("filter");
        final List<String> filterArgs = new ArrayList<>();
        filterArgs.add(prefix+".fa");
        final SGAModule.RuntimeInfo filterInfo = filter.run(sgaPath, outputDirectory, filterArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(filterInfo);

        return prefix+".filter.pass.fa";
    }

    @VisibleForTesting
    static String runSGARmDuplicate(final Path sgaPath,
                                    final File inputFASTAFile,
                                    final File outputDirectory,
                                    final SGAModule indexer,
                                    final List<String> indexerArgs,
                                    final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                    final boolean enableSTDIOCapture){
        return runSimpleModuleFollowedByIndexing(sgaPath, "rmdup", ".rmdup.fa", inputFASTAFile, outputDirectory, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
    }

    @VisibleForTesting
    static String runSGAFMMerge(final Path sgaPath,
                                final File inputFASTAFile,
                                final File outputDirectory,
                                final SGAModule indexer,
                                final List<String> indexerArgs,
                                final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                final boolean enableSTDIOCapture){
        return runSimpleModuleFollowedByIndexing(sgaPath, "fm-merge", ".merged.fa", inputFASTAFile, outputDirectory, indexer, indexerArgs, collectedRuntimeInfo, enableSTDIOCapture);
    }

    @VisibleForTesting
    static String runSGAOverlapAndAssemble(final Path sgaPath,
                                           final File inputFASTAFile,
                                           final File outputDirectory,
                                           final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                           final boolean enableSTDIOCapture){

        final SGAModule overlap = new SGAModule("overlap");
        final List<String> overlapArgs = new ArrayList<>();
        overlapArgs.add(inputFASTAFile.getName());

        final SGAModule.RuntimeInfo overlapInfo = overlap.run(sgaPath, outputDirectory, overlapArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(overlapInfo);

        final String prefix = FilenameUtils.getBaseName(inputFASTAFile.getName());

        final SGAModule assemble = new SGAModule("assemble");
        final List<String> assembleArgs = new ArrayList<>();
        assembleArgs.add("--out-prefix"); assembleArgs.add(prefix);
        assembleArgs.add(prefix+".asqg.gz");
        final SGAModule.RuntimeInfo assembleInfo = assemble.run(sgaPath, outputDirectory, assembleArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(assembleInfo);

        return prefix+"-contigs.fa";
    }

    // boiler plate code for running simple sga modules (simple in the sense that no options needs to be specified to make it work)
    // that perform a task followed by indexing its output
    private static String runSimpleModuleFollowedByIndexing(final Path sgaPath,
                                                            final String moduleName,
                                                            final String extensionToAppend,
                                                            final File inputFASTAFile,
                                                            final File outputDirectory,
                                                            final SGAModule indexer,
                                                            final List<String> indexerArgs,
                                                            final List<SGAModule.RuntimeInfo> collectedRuntimeInfo,
                                                            final boolean enableSTDIOCapture){

        final SGAModule module = new SGAModule(moduleName);
        final List<String> args = new ArrayList<>();
        args.add(inputFASTAFile.getName());

        final SGAModule.RuntimeInfo moduleInfo = module.run(sgaPath, outputDirectory, args, enableSTDIOCapture);
        collectedRuntimeInfo.add(moduleInfo);

        final String outputFileName = FilenameUtils.getBaseName(inputFASTAFile.getName()) + extensionToAppend;

        indexerArgs.set(indexerArgs.size()-1, outputFileName);
        final SGAModule.RuntimeInfo indexerInfo = indexer.run(sgaPath, outputDirectory, indexerArgs, enableSTDIOCapture);
        collectedRuntimeInfo.add(indexerInfo);

        return outputFileName;
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

    /**
     * Represents a collection of assembled contigs (not including the variants) produced by "sga assemble".
     */
    @VisibleForTesting
    static final class ContigsCollection implements Serializable{
        private static final long serialVersionUID = 1L;

        @VisibleForTesting
        static final class ContigSequence implements Serializable{
            private static final long serialVersionUID = 1L;

            private final String sequence;
            public ContigSequence(final String sequence){ this.sequence = sequence; }

            @Override
            public String toString(){
                return sequence;
            }
        }

        @VisibleForTesting
        static final class ContigID implements Serializable{
            private static final long serialVersionUID = 1L;

            private final String id;
            public ContigID(final String idString) { this.id = idString; }

            @Override
            public String toString(){
                return id;
            }
        }

        private final List<Tuple2<ContigID, ContigSequence>> contents;

        public List<Tuple2<ContigID, ContigSequence>> getContents(){
            return contents;
        }

        @Override
        public String toString(){
            return StringUtils.join(toListOfStrings(),"\n");
        }

        public String toOneLineString(){
            return StringUtils.join(toListOfStrings(),"|");
        }

        public List<String> toListOfStrings(){
            if(null==contents){
                return null;
            }
            final List<String> res = new ArrayList<>();
            for(final Tuple2<ContigID, ContigSequence> contig : contents){
                res.add(contig._1().toString());
                res.add(contig._2().toString());
            }
            return res;
        }

        public ContigsCollection(final List<String> fileContents){

            if(null==fileContents){
                contents = null;
            }else{
                contents = new ArrayList<>();
                for(int i=0; i<fileContents.size(); i+=2){
                    contents.add(new Tuple2<>(new ContigID(fileContents.get(i)), new ContigSequence(fileContents.get(i+1))));
                }
            }
        }
    }
}