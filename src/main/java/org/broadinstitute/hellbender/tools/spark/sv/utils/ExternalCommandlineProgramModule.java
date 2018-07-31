package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.base.Throwables;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * A class for experimenting with cmd line programs (such as bwa) that's not part of GATK (or has Java bindings yet).
 * See {@code ExternalCommandlineProgramModule.SGAModule} for illustration on how to construct a module,
 * and {@code ExternalCommandlineProgramModule#exampleRun(String, String, File, boolean)}] on how to execute it.
 */
abstract class ExternalCommandlineProgramModule {

    public ExternalCommandlineProgramModule(){  }

    /**
     * Runtime information of this module collected throughout its run.
     */
    public static final class RuntimeInfo implements Serializable{
        private static final long serialVersionUID = 1L;

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
        public final ReturnStatus returnStatus;
        public final String stdoutMsg;
        public final String stderrMsg;
        public final String moduleName;

        public RuntimeInfo(final String moduleName,
                           final ReturnStatus status,
                           final String stdoutMsg,
                           final String stderrMsg){

            this.moduleName   = moduleName;
            this.returnStatus = status;
            this.stdoutMsg    = stdoutMsg;
            this.stderrMsg    = stderrMsg;
        }

        // for parsing contents to human readable string
        @Override
        public String toString(){
            String result = moduleName + "\t";
            switch (returnStatus){
                case STARTFAIL:
                    result += " failed to start.\n"; break;
                case INTERRUPTION:
                    result += " was interrupted.\n"; break;
                case STDIOFAIL:
                    result += " stdout and stderr wasn't successfully captured.\n"; break;
                case PGFAIL:
                    result += " returned with non-zero exit status, see stderr message for detailed return status.\n"; break;
                default://SUCCESS
                    result += " successfully executed.\n";
            }
            result += " Module stdout message: \n";
            result += null!=stdoutMsg ? stdoutMsg+"\n" : "is empty\n";

            result += " Module stderr message: \n";
            result += null!=stderrMsg ? stderrMsg+"\n" : "is empty\n";
            return result;
        }
    }

    abstract public String getModuleName();

    /**
     * Builds the initial argument list for the command line (e.g. "bwa mem").
     * @return a list of strings having the names of requested tools
     */
    abstract public List<String> initializeCommands(final Path pathToProgram);

    /**
     * Function to execute the command line program, if available, on the provided arguments.
     * @param pathToProgram             full path to the underlying program
     * @param directoryToWorkIn         a directory for the program to work in (mainly for output); can be null
     * @param runtimeArguments          arguments to be provided to the program; can be null (e.g. program like "ls" taking no arguments)
     * @param workingEnvironmentArgs    arguments to setup working environment for the program to run under
     * @return                          exit status of running the process, paired with stdout and stderr strings of the process
     *                                  except for case SUCCESS, accompanying stdout is null, and stderr contains error message
     */
    public RuntimeInfo run(final Path pathToProgram,
                           final File directoryToWorkIn,
                           final List<String> runtimeArguments,
                           final boolean enableSTDIOCapture,
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
            final File stdoutFile = enableSTDIOCapture ? Files.createTempFile("stdout", "log").toFile() : null;
            final File stderrFile = enableSTDIOCapture ? Files.createTempFile("stderr", "log").toFile() : null;
            if(enableSTDIOCapture){
                builder.redirectOutput(stdoutFile);
                builder.redirectError(stderrFile);
            }

            String out = "";
            String err = "";

            final Process process = builder.start();
            final int exitStatus = process.waitFor();

            if(enableSTDIOCapture){
                try{
                    out = FileUtils.readFileToString(stdoutFile, Charset.defaultCharset());
                    err = FileUtils.readFileToString(stderrFile, Charset.defaultCharset());
                }catch (final IOException ex){
                    return new RuntimeInfo(getModuleName(), RuntimeInfo.ReturnStatus.STDIOFAIL, null, ex.getMessage());
                }
            }

            if(0!=exitStatus){
                return new RuntimeInfo(getModuleName(),
                                       RuntimeInfo.ReturnStatus.PGFAIL,
                                       null,
                                       "Program returned exit status " + exitStatus + "\n" + String.join(" ", commands) + "\n" + out + "\n" + err);
            }

            return new RuntimeInfo(getModuleName(), RuntimeInfo.ReturnStatus.SUCCESS, out, err);
        } catch (final IOException e){
            return new RuntimeInfo(getModuleName(), RuntimeInfo.ReturnStatus.STARTFAIL, null, e.getMessage() + Throwables.getStackTraceAsString(e));
        } catch (final InterruptedException e){
            return new RuntimeInfo(getModuleName(), RuntimeInfo.ReturnStatus.INTERRUPTION, null, e.getMessage() + Throwables.getStackTraceAsString(e));
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

    /**
     * INTENDED FOR ILLUSTRATION PURPOSE ONLY FOR HOW TO CONSTRUCT A {@link ExternalCommandlineProgramModule}
     * Represents an SGA module that can be called via "run" (in base) to do actual work.
     */
    @Deprecated
    private static final class SGAModule extends ExternalCommandlineProgramModule {

        private final String moduleName;

        SGAModule(final String moduleName){
            this.moduleName = moduleName;
        }

        @Override
        public String getModuleName(){
            return "sga " + moduleName;
        }

        @Override
        public List<String> initializeCommands(final Path pathToSGA) {
            final ArrayList<String> result = new ArrayList<>();
            result.add(pathToSGA.toString());
            result.add(moduleName);
            return result;
        }
    }

    /**
     * INTENDED FOR ILLUSTRATION PURPOSE ONLY FOR HOW TO USE A {@link ExternalCommandlineProgramModule}
     * @param sgaPathString         absolute path to the program
     * @param moduleName            name of the module to be run (e.g. "bwa index", "bwa mem")
     * @param workingDir            directory to work in, especially if there's file IO
     * @param enableSTDIOCapture    whether to enable capturing stdout and stderr from the module
     */
    @Deprecated
    private static void exampleRun(final String sgaPathString,
                                   final String moduleName,
                                   final File workingDir,
                                   final boolean enableSTDIOCapture) {
        // intake program path
        final Path sgaPath = Paths.get(sgaPathString);
        // instantiate the desired module, and provide appropriate parameters
        final SGAModule module = new SGAModule(moduleName);
        final List<String> moduleArgs = new ArrayList<>();
        moduleArgs.add("--algorithm"); moduleArgs.add("ropebwt");
        moduleArgs.add("--check");
        moduleArgs.add("");
        // run, collect runtime information and check return status
        final SGAModule.RuntimeInfo moduleRuntimeInfo = module.run(sgaPath, workingDir, moduleArgs, enableSTDIOCapture);
        final SGAModule.RuntimeInfo.ReturnStatus returnStatus = moduleRuntimeInfo.returnStatus;
        if(!(returnStatus.equals(SGAModule.RuntimeInfo.ReturnStatus.SUCCESS))){
            throw new RuntimeException("Appropriate action and error message here");
        }
    }
}
