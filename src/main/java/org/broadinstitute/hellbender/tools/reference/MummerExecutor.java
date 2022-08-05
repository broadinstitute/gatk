package org.broadinstitute.hellbender.tools.reference;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.File;
import java.nio.file.Path;
import java.util.*;


public final class MummerExecutor {

    private static final Logger logger = LogManager.getLogger(MummerExecutor.class);
    private File mummerExecutableDirectory;

    public MummerExecutor(File mummerExecutableDirectory){
        this.mummerExecutableDirectory = mummerExecutableDirectory;
    }

    public File getMummerExecutableDirectory() {
        return mummerExecutableDirectory;
    }

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
        File showSNPSOutput = IOUtils.createTempFile("showSNPSOutput", ".snps");
        String[] showSNPsArgs = {mummerExecutableDirectory.getAbsolutePath() + "/show-snps", "-rlTH", deltaFilterOutput.getAbsolutePath()};
        ProcessOutput showSNPs = runShellCommand(showSNPsArgs, null, showSNPSOutput, false);

        // ALL2VCF
        logger.debug("Running all2vcf.");
        //File tempVCF = IOUtils.createTempFileInDirectory("tempVCF", ".vcf", outputDirectory);
        File tempVCF = new File(outputDirectory, "testVCF.vcf");
        String script = "/Users/ocohen/workingcode/MUMmer3.23/all2vcf";
        List<String> all2vcfArgs = Arrays.asList("--snps", showSNPSOutput.getAbsolutePath(), "--reference", fasta1.getAbsolutePath(), "--output-header");
        ProcessOutput all2vcf = runPythonCommand(script, all2vcfArgs, null, tempVCF, false);

        return tempVCF;
    }

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

}
