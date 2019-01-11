package org.broadinstitute.hellbender.utils.R;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.runtime.ScriptExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Generic service for executing RScripts
 */
public final class RScriptExecutor extends ScriptExecutor {
    private static final Logger logger = LogManager.getLogger(RScriptExecutor.class);
    private final List<RScriptLibrary> libraries = new ArrayList<>();
    private final List<Resource> scriptResources = new ArrayList<>();
    private final List<File> scriptFiles = new ArrayList<>();
    private final List<String> args = new ArrayList<>();

    public RScriptExecutor() {
        this(false);
    }

    /**
     * @param ensureExecutableExists throw if the RScript executable cannot be located
     */
    public RScriptExecutor(boolean ensureExecutableExists) {
        super("Rscript");
        if (ensureExecutableExists && !externalExecutableExists()) {
            executableMissing();
        }
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

    @Override
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

    /**
     * Return an exception specific to this executor type, to be thrown on error conditions.
     * @param message
     */
    @Override
    public RScriptExecutorException getScriptException(final String message) {
        return new RScriptExecutorException(message.toString());
    }

    public boolean exec() {
        final List<File> tempDirs = new ArrayList<>();
        try {
            File tempLibSourceDir  = IOUtils.createTempDir("RlibSources.");
            File tempLibInstallationDir = IOUtils.createTempDir("Rlib.");
            tempDirs.add(tempLibSourceDir);
            tempDirs.add(tempLibInstallationDir);

            StringBuilder expression = new StringBuilder("tempLibDir = '").append(tempLibInstallationDir).append("';");

            if (!this.libraries.isEmpty()) {
                List<String> tempLibraryPaths = new ArrayList<>();
                for (RScriptLibrary library: this.libraries) {
                    final File tempLibrary = library.writeLibraryToTempFile(tempLibSourceDir);
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
                final File tempScript = IOUtils.writeTempResource(script);
                expression.append("source('").append(tempScript.getAbsolutePath()).append("');");
            }

            for (File script: this.scriptFiles) {
                expression.append("source('").append(script.getAbsolutePath()).append("');");
            }

            String[] cmd = new String[this.args.size() + 3];
            int i = 0;
            cmd[i++] = externalScriptExecutableName;
            cmd[i++] = "-e";
            cmd[i++] = expression.toString();
            for (String arg: this.args)
                cmd[i++] = arg;

            // actually run the script
            return executeCuratedArgs(cmd);
        } catch (GATKException e) {
            if (!ignoreExceptions) {
                throw e;
            } else {
                logger.warn(e.getMessage());
                return false;
            }
        } finally {
            for (final File tempDir: tempDirs) {
                // deletes the dirs and their contents
                FileUtils.deleteQuietly(tempDir);
            }
        }
    }
}
