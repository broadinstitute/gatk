/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.runtime.RuntimeUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Generic service for executing RScripts
 */
public class RScriptExecutor {
    private static final String RSCRIPT_BINARY = "Rscript";
    private static final File RSCRIPT_PATH = RuntimeUtils.which(RSCRIPT_BINARY);
    public static final boolean RSCRIPT_EXISTS = (RSCRIPT_PATH != null);
    private static final String RSCRIPT_MISSING_MESSAGE = "Please add the Rscript directory to your environment ${PATH}";

    /**
     * our log
     */
    private static Logger logger = LogManager.getLogger(RScriptExecutor.class);

    private boolean exceptOnError = false;
    private final List<RScriptLibrary> libraries = new ArrayList<RScriptLibrary>();
    private final List<Resource> scriptResources = new ArrayList<Resource>();
    private final List<File> scriptFiles = new ArrayList<File>();
    private final List<String> args = new ArrayList<String>();

    public void setExceptOnError(boolean exceptOnError) {
        this.exceptOnError = exceptOnError;
    }

    public void addLibrary(RScriptLibrary library) {
        this.libraries.add(library);
    }

    public void addScript(Resource script) {
        this.scriptResources.add(script);
    }

    public void addScript(File script) {
        this.scriptFiles.add(script);
    }

    /**
     * Adds args to the end of the Rscript command line.
     * @param args the args.
     * @throws NullPointerException if any of the args are null.
     */
    public void addArgs(Object... args) {
        for (Object arg: args)
            this.args.add(arg.toString());
    }

    public String getApproximateCommandLine() {
        StringBuilder command = new StringBuilder("Rscript");
        for (Resource script: this.scriptResources)
            command.append(" (resource)").append(script.getFullPath());
        for (File script: this.scriptFiles)
            command.append(" ").append(script.getAbsolutePath());
        for (String arg: this.args)
            command.append(" ").append(arg);
        return command.toString();
    }

    public boolean exec() {
        if (!RSCRIPT_EXISTS) {
            if (exceptOnError) {
                throw new UserException.CannotExecuteRScript(RSCRIPT_MISSING_MESSAGE);
            } else {
                logger.warn("Skipping: " + getApproximateCommandLine());
                return false;
            }
        }

        List<File> tempFiles = new ArrayList<>();
        try {
            File tempLibSourceDir  = IOUtils.tempDir("RlibSources.", "");
            File tempLibInstallationDir = IOUtils.tempDir("Rlib.", "");
            tempFiles.add(tempLibSourceDir);
            tempFiles.add(tempLibInstallationDir);

            StringBuilder expression = new StringBuilder("tempLibDir = '").append(tempLibInstallationDir).append("';");

            if (this.libraries.size() > 0) {
                List<String> tempLibraryPaths = new ArrayList<String>();
                for (RScriptLibrary library: this.libraries) {
                    File tempLibrary = library.writeLibrary(tempLibSourceDir);
                    tempFiles.add(tempLibrary);
                    tempLibraryPaths.add(tempLibrary.getAbsolutePath());
                }

                expression.append("install.packages(");
                expression.append("pkgs=c('").append(StringUtils.join(tempLibraryPaths, "', '")).append("'), lib=tempLibDir, repos=NULL, type='source', ");
                // Install faster by eliminating cruft.
                expression.append("INSTALL_opts=c('--no-libs', '--no-data', '--no-help', '--no-demo', '--no-exec')");
                expression.append(");");

                for (RScriptLibrary library: this.libraries) {
                    expression.append("library('").append(library.getLibraryName()).append("', lib.loc=tempLibDir);");
                }
            }

            for (Resource script: this.scriptResources) {
                File tempScript = IOUtils.writeTempResource(script);
                tempFiles.add(tempScript);
                expression.append("source('").append(tempScript.getAbsolutePath()).append("');");
            }

            for (File script: this.scriptFiles) {
                expression.append("source('").append(script.getAbsolutePath()).append("');");
            }

            String[] cmd = new String[this.args.size() + 3];
            int i = 0;
            cmd[i++] = RSCRIPT_BINARY;
            cmd[i++] = "-e";
            cmd[i++] = expression.toString();
            for (String arg: this.args)
                cmd[i++] = arg;

            ProcessSettings processSettings = new ProcessSettings(cmd);
            if (logger.isDebugEnabled()) {
                processSettings.getStdoutSettings().printStandard(true);
                processSettings.getStderrSettings().printStandard(true);
            }

            ProcessController controller = ProcessController.getThreadLocal();

            if (logger.isDebugEnabled()) {
                logger.debug("Executing:");
                for (String arg: cmd)
                    logger.debug("  " + arg);
            }
            int exitValue = controller.exec(processSettings).getExitValue();
            logger.debug("Result: " + exitValue);

            if (exitValue != 0)
                throw new RScriptExecutorException(
                        "RScript exited with " + exitValue +
                                (logger.isDebugEnabled() ? "" : ". Run with -l DEBUG for more info."));

            return true;
        } catch (GATKException e) {
            if (exceptOnError) {
                throw e;
            } else {
                logger.warn(e.getMessage());
                return false;
            }
        } finally {
            for (File temp: tempFiles)
                FileUtils.deleteQuietly(temp);
        }
    }
}
