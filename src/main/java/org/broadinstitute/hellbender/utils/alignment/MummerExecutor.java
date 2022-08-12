package org.broadinstitute.hellbender.utils.alignment;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.SystemUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.NativeUtils;
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
 * Class for executing MUMmer alignment pipeline to detect SNPs and INDELs in mismatching sequences.
 * The MUMmer pipeline runs nucmer, delta-filter, and show-snps programs, sequentially.
 * Used in CompareReferences tool when running in FULL_ALIGNMENT mode.
 */
public final class MummerExecutor {

    public static final String MUMMER_X86_64_MAC_BINARIES_ZIPFILE = "MUMmer-4.0.0rc1_mac_x64_64.zip";
    public static final String MUMMER_X86_64_LINUX_BINARIES_ZIPFILE = "MUMmer-4.0.0rc1_linux_x64_64.zip";
    // TODO: Need to recompile MUMmer 4.0.0rc1 for native M1 architecture as well
    // TODO: (however the Mac Intel build should be runnable on M1 platforms with x86 emulation)

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
        final Resource mummerDistribution = new Resource(selectCorrectMummerDistributionForPlatform(), getClass());
        mummerExecutableDirectory = prepareMUMmerExecutionDirectory(mummerDistribution);
    }

    private String selectCorrectMummerDistributionForPlatform() {
        if ( NativeUtils.runningOnMac() && NativeUtils.runningOn64BitX86Architecture() ) {
            logger.info("Using the x86_64 Mac OS build of MUMmer");
            return MUMMER_X86_64_MAC_BINARIES_ZIPFILE;
        }
        else if ( NativeUtils.runningOnLinux() && NativeUtils.runningOn64BitX86Architecture() ) {
            logger.info("Using the x86_64 Linux build of MUMmer");
            return MUMMER_X86_64_LINUX_BINARIES_ZIPFILE;
        }
        else {
            throw new UserException("Unable to run MUMmer aligner: GATK does not contain a MUMmer distribution for system architecture " +
                    SystemUtils.OS_ARCH + " on operating system " + SystemUtils.OS_NAME);
        }
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
    public File executeMummer(File fasta1, File fasta2, File outputDirectory){

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
        File showSNPSOutput = new File(outputDirectory, "snps_output.snps");
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

    // method to unzip and locate the MUMmer binaries packaged within GATK
    @VisibleForTesting
    static File prepareMUMmerExecutionDirectory(final Resource mummerZipFile){
        try{
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
