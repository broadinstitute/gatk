package org.broadinstitute.hellbender.utils.alignment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.PosixFilePermissions;
import java.util.*;

/**
 * Class for executing MUMmer pipeline.
 */
public final class MummerExecutor {

    // this is an M1 specific build of MUMmer
    // TODO: add intel build
    public static final String MUMMER_BINARIES_ZIPFILE = "MUMmer3.23Binaries.zip";
    private static final Logger logger = LogManager.getLogger(MummerExecutor.class);
    private File mummerExecutableDirectory;

    /**
     * Returns a MummerExecutor pointing to the provided directory
     * @param mummerExecutableDirectory File representing the location of the MUMmer executables
     */
    public MummerExecutor(File mummerExecutableDirectory){
        this.mummerExecutableDirectory = mummerExecutableDirectory;
    }


    /**
     * Returns a MummerExecutor pointing to the unzipped MUMmer executables packaged in GATK
     */
    public MummerExecutor(){
        mummerExecutableDirectory = prepareMUMmerExecutionDirectory();
    }

    /**
     * Getter method to access the location of the MUMmer executables.
     * @return File representing the directory location of the MUMmer executables
     */
    public File getMummerExecutableDirectory() {
        return mummerExecutableDirectory;
    }

    /**
     * Execute a full alignment of a sequence across a given pair of references to get a File containing the SNPs and INDELs.
     * Begins by running nucmer to get a delta file, which is passed to delta-filter. The result of delta filter is a
     * delta file filtered to only output alignments of length >= 1; this file is passed to show-snps, which displays SNPs
     * and INDELs in a tabular format.
     *
     * @param fasta1 first reference
     * @param fasta2 second reference
     * @param outputDirectory directory to output final snps file
     * @return the final snps File
     */
    public File executeMummer(File fasta1, File fasta2, File outputDirectory, String sequenceName){

        // NUCMER
        logger.debug("Running nucmer.");
        File nucmerTempDirectory = IOUtils.createTempDir("nucmerTempDir");
        File deltaFile = new File(nucmerTempDirectory, "deltaFile"); // delta file for nucmer output --> input to delta-filter
        String[] nucmerArgs = {mummerExecutableDirectory.getAbsolutePath() + "/nucmer", "--mum", "-p", deltaFile.getAbsolutePath(), fasta1.getAbsolutePath(), fasta2.getAbsolutePath()};
        ProcessOutput nucmer = runShellCommand(nucmerArgs, null, null,false);

        // DELTA-FILTER
        logger.debug("Running delta-filter.");
        File deltaFilterOutput = IOUtils.createTempFile("deltaFilterOutput", ".delta"); // file for delta filter output --> input to show-snps
        String[] deltaFilterArgs = {mummerExecutableDirectory.getAbsolutePath() + "/delta-filter", "-1", deltaFile.getAbsolutePath() + ".delta"};
        ProcessOutput deltaFilter = runShellCommand(deltaFilterArgs, null, deltaFilterOutput, false);

        // SHOW-SNPS
        logger.debug("Running show-snps.");
        File showSNPSOutput = new File(outputDirectory, String.format("chr%s_snps_output.snps", sequenceName));
        String[] showSNPsArgs = {mummerExecutableDirectory.getAbsolutePath() + "/show-snps", "-rlTH", deltaFilterOutput.getAbsolutePath()};
        ProcessOutput showSNPs = runShellCommand(showSNPsArgs, null, showSNPSOutput, false);

        return showSNPSOutput;
    }

    /**
     * Method to run shell commands from GATK.
     *
     * @param command String[] containing the commandline for execution
     * @param environment a Map of environment variables to their values; may be left null if no environment variables needed
     * @param stdoutCaptureFile File for output
     * @param printStdout boolean to trigger displaying to stdout
     * @return the ProcessOutput of the run
     */
    public static ProcessOutput runShellCommand(String[] command, Map<String, String> environment, File stdoutCaptureFile, boolean printStdout){
        ProcessController processController = ProcessController.getThreadLocal();
        final ProcessSettings prs = new ProcessSettings(command);
        if(printStdout){
            prs.getStderrSettings().printStandard(true);
            prs.getStdoutSettings().printStandard(true);
        }
        if(stdoutCaptureFile != null){
            prs.getStdoutSettings().setOutputFile(stdoutCaptureFile);
        }
        prs.setEnvironment(environment);
        final ProcessOutput output = processController.exec(prs);

        if(output.getExitValue() != 0){
            throw new UserException("Error running " + command[0] + ". Exit with code " + output.getExitValue());
        }

        return output;
    }

    /**
     * Method to run python commands from GATK.
     *
     * @param script the name of the executable as a String
     * @param scriptArguments a List of Strings containing the remainder of the commandline for execution
     * @param additionalEnvironmentVars a Map of environment variables to their values; may be left null if no environment variables needed
     * @param stdoutCaptureFile File for output
     * @param printStdout boolean to trigger displaying to stdout
     * @return the ProcessOutput of the run
     */
    public static ProcessOutput runPythonCommand(String script, List<String> scriptArguments, Map<String, String> additionalEnvironmentVars, File stdoutCaptureFile, boolean printStdout){
        Map<String, String> environment = new HashMap<>();
        environment.putAll(System.getenv());
        if(additionalEnvironmentVars != null){
            environment.putAll(additionalEnvironmentVars);
        }
        List<String> args = new ArrayList<>();
        args.add("python");
        args.add(script);
        args.addAll(scriptArguments);
        ProcessOutput output = runShellCommand(args.toArray(new String[]{}), environment, stdoutCaptureFile, printStdout);
        if(output.getExitValue() != 0){
            throw new UserException("Error running " + script + ". Exit with code " + output.getExitValue());
        }
        return output;
    }

    // method to locate the MUMmer binaries packaged within GATK
    private File prepareMUMmerExecutionDirectory(){
        try{
            Resource mummerZipFile = new Resource(MUMMER_BINARIES_ZIPFILE, getClass());
            File tempMUMmerZipFile = IOUtils.writeTempResource(mummerZipFile);
            File mummerExecutionDirectory = IOUtils.createTempDir("MUMmerExecutionDirectory");
            IOUtils.unzipToFolder(Paths.get(tempMUMmerZipFile.getAbsolutePath()), Paths.get(mummerExecutionDirectory.getAbsolutePath()));
            for(File file : mummerExecutionDirectory.listFiles()){
                Path path = Paths.get(file.getAbsolutePath());
                Files.setPosixFilePermissions(path, PosixFilePermissions.fromString("rwxr-xr-x"));
            }
            return mummerExecutionDirectory;
        }
        catch(IOException e){
            throw new UserException("Unable to unzip MUMmer binaries.", e);
        }
    }

}
