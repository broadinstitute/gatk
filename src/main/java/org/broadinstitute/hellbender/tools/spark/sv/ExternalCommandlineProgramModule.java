package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.base.Throwables;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import scala.Tuple3;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * A class for experimenting with cmd line programs (such as bwa and sga) that's not part of GATK (or has Java bindings yet).
 */
public abstract class ExternalCommandlineProgramModule {

    public ExternalCommandlineProgramModule(){  }

    /**
     * Return status of running the underlying program:
     * STARTFAIL:       program failed to start
     * INTERRUPTION:    program was interrupted
     * STDIOFAIL:       failed to capture stdout and stderr of the program
     * PGRAIL:          program failed for some other reasons; error messages will be logged separately
     */
    public enum ReturnStatus{
        SUCCESS, STARTFAIL, INTERRUPTION, STDIOFAIL, PGFAIL
    }

    /**
     * Builds the initial argument list for the command line (e.g. "bwa mem").
     * @return a list of strings having the names of requested tools
     */
    abstract List<String> initializeCommands(final Path pathToProgram);

    /**
     * Function to execute the command line program, if available, on the provided arguments.
     * @param pathToProgram             full path to the underlying program
     * @param directoryToWorkIn         a directory for the program to work in (mainly for output); can be null
     * @param runtimeArguments          arguments to be provided to the program; can be null (e.g. program like "ls" taking no arguments)
     * @param workingEnvironmentArgs    arguments to setup working environment for the program to run under
     * @return                          exit status of running the process, paired with stdout and stderr strings of the process
     *                                  except for case SUCCESS, accompanying stdout is null, and stderr contains error message
     */
    public Tuple3<ReturnStatus, String, String> run(final Path pathToProgram,
                                                    final File directoryToWorkIn,
                                                    final List<String> runtimeArguments,
                                                    final String... workingEnvironmentArgs) {

        List<String> commands = initializeCommands(pathToProgram);
        if(null!=runtimeArguments && !runtimeArguments.isEmpty()){
            commands.addAll(runtimeArguments);
        }

        final ProcessBuilder builder = new ProcessBuilder(commands);

        setupWorkingEnvironment(builder, workingEnvironmentArgs);
        if(null!=directoryToWorkIn){
            builder.directory(directoryToWorkIn);
        }

        try{
            final File stdoutFile = Files.createTempFile("stdout", "log").toFile();
            final File stderrFile = Files.createTempFile("stderr", "log").toFile();
            builder.redirectOutput(stdoutFile);
            builder.redirectError(stderrFile);

            String out = "";
            String err = "";

            final Process process = builder.start();
            final int exitStatus = process.waitFor();

            try{
                out = FileUtils.readFileToString(stdoutFile, Charset.defaultCharset());
                err = FileUtils.readFileToString(stderrFile, Charset.defaultCharset());
            }catch (final IOException ex){
                return new Tuple3<>(ReturnStatus.STDIOFAIL, null, ex.getMessage());
            }

            if(0!=exitStatus){
                return new Tuple3<>(ReturnStatus.PGFAIL,
                                    null,
                                    "Program returned exit status " + exitStatus + "\n" + String.join(" ", commands) + "\n" + out + "\n" + err);
            }

            return new Tuple3<>(ReturnStatus.SUCCESS, out, err);
        } catch (final IOException e){
            return new Tuple3<>(ReturnStatus.STARTFAIL, null, e.getMessage() + Throwables.getStackTraceAsString(e));
        } catch (final InterruptedException e){
            return new Tuple3<>(ReturnStatus.INTERRUPTION, null, e.getMessage() + Throwables.getStackTraceAsString(e));
        }
    }

    /**
     * Sets up working environment for the program.
     * Inheriting classes can, and most likely should, override to do non-trivial work.
     * @param builder process builder for the program
     * @param args    arguments used for setting up the environment
     */
    protected void setupWorkingEnvironment(final ProcessBuilder builder, final String... args){

    }
}
