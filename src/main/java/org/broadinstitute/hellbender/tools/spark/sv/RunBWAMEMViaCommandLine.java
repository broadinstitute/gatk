package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import scala.Tuple3;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.ArrayList;

// TODO: choose which parameters allowed to be tunable
// TODO: if throws, would temp files be cleaned up automatically?
@CommandLineProgramProperties(
        summary        = "Minimal program to call BWAMEM for performing alignment, allowing limited options.",
        oneLineSummary = "Minimal work to call BWAMEM.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class RunBWAMEMViaCommandLine extends CommandLineProgram {

    @Argument(doc       = "Absolute path to BWA program.",
              shortName = "bwaPath",
              fullName  = "fullPathToBWA",
              optional  = false)
    public String pathToBWA = null;

    @PositionalArguments(minElements = 1, maxElements = 2, doc = "Path to input FASTQ/A file(s) to be aligned.")
    public List<File> inputFastaFiles;

    @Argument(doc       = "If set to true, performs smart pairing (input FASTQ assumed interleaved), and ignores second fastq input.",
              shortName = "p",
              fullName  = "interLeaved",
              optional  = true)
    public boolean interLeaved = false;

    @Argument(doc       = "If set to true, assumes input are single ended data.",
              shortName = "se",
              fullName  = "singleEnd",
              optional  = true)
    public boolean SEInput = false;

    @Argument(doc       = "Path to reference of the target organism, if alignment of assembled contigs to reference is desired.",
              shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
              fullName  = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
              optional  = false)
    public String reference = null;

    @Argument(doc       = "Sam file to write results to.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              optional  = false)
    public String samOutput = null;

    @Argument(doc       = "File name where stderr of running bwa should be redirected to.",
              shortName = "eFile",
              fullName  = "stderFile",
              optional  = true)
    public String stderrDestFileName = null;

    @Argument(doc       = "Number of threads to use when running bwa.",
              shortName = "t",
              fullName  = "threads",
              optional  = true)
    public int threads = 1;

    @Argument(doc       = "Number of threads to use when running bwa.",
              shortName = "K",
              fullName  = "chunkSizeEachThread",
              optional  = true)
    public long chunkSize = 0;

    /**
     * Validates user options. Throws UserException if arguments provided don't make logical sense.
     * Runs bwa mem. Throws GATKException if the bwa mem process erred.
     * @return  stderr message from underlying program (e.g. logs, performance, etc) upon successful execution,
     *          empty if caller decides to redirect the stderr message to a file
     */
    @Override
    public String doWork(){

        validateUserOptions();

        final BWAMEMModule bwamem = new BWAMEMModule();
        final Tuple3<ExternalCommandlineProgramModule.ReturnStatus, String, String> result = bwamem.run(Paths.get(pathToBWA),
                                                                                                        new File(System.getProperty("user.dir")),
                                                                                                        makeArgs());

        return validateResults(result);
    }

    private void validateUserOptions() throws UserException {

        inputFastaFiles.forEach(IOUtil::assertFileIsReadable);

        final boolean validOptions = (inputFastaFiles.size()==1) ? !(interLeaved && SEInput) : (!interLeaved || SEInput);
        if(!validOptions){
            throw new UserException("CMD line argument options on input (paired, interleaved, or SE) don't make sense. Please check.");
        }
    }

    private String validateResults(final Tuple3<ExternalCommandlineProgramModule.ReturnStatus, String, String> result) throws GATKException {

        switch (result._1()){
            case STARTFAIL:     throw new GATKException("Failed to start bwa mem.\n" + result._3());
            case INTERRUPTION:  throw new GATKException("The bwa mem process was interrupted.\n" + result._3());
            case STDIOFAIL:     throw new GATKException("Failed to capture bwa stdout/stderr message\n" + result._3());
            case PGFAIL:        throw new GATKException(result._3());
            default:            return writeSamFile(result);
        }
    }

    /**
     * Writes stdout from bwa mem to designated file.
     * Return stderr message as string, which will be empty if user decides to redirect stderr message to file (i.e. stderrDestFileName is set)
     */
    private String writeSamFile(final Tuple3<ExternalCommandlineProgramModule.ReturnStatus, String, String> result){

        try{
            final File samFile = new File(samOutput);
            samFile.createNewFile();
            FileUtils.writeStringToFile(samFile, result._2());

        } catch (final IOException e){
            throw new GATKException("Failed to dump bwa mem stdout to designated sam file.\n"  + e.getMessage());
        }

        if(null==stderrDestFileName){
            return result._3();
        }else{
            final File stderrFile = new File(stderrDestFileName);
            try{
                FileUtils.writeStringToFile(stderrFile, result._3());
            }catch(final IOException e){
                System.err.println("Failed to dump stderr message form bwa mem to designated file:\n" + stderrDestFileName);
            }
            return result._3();
        }
    }

    private List<String> makeArgs(){
        final List<String> args = new ArrayList<>();

        final Path pathToReference = Paths.get(reference);
        final Path pathToInput = inputFastaFiles.get(0).toPath();

        args.add("-t");
        args.add(Integer.toString(threads));

        if(0!=chunkSize){
            args.add("-K");
            args.add(Long.toString(this.chunkSize));
        }

        if(interLeaved){ // paired reads, interleaved
            args.add("-p");
        }else if(1==inputFastaFiles.size()) { // SE reads
            args.add("-S"); // skips mate rescuing and pairing
            args.add("-P");
        }
        args.add(pathToReference.toAbsolutePath().toString());
        args.add(pathToInput.toAbsolutePath().toString());

        if(2==inputFastaFiles.size()){ // paired reads, separate files
            final Path pathToSecondInput = inputFastaFiles.get(1).toPath().toAbsolutePath();
            args.add(pathToSecondInput.toString());
        }

        return args;
    }
}