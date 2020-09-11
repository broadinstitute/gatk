package org.broadinstitute.hellbender.utils.runtime;

public final class ProcessOutput {
    private final int exitValue;
    private final StreamOutput stdout;
    private final StreamOutput stderr;

    /**
     * The output of a process.
     *
     * @param exitValue The exit value.
     * @param stdout    The capture of stdout as defined by the stdout OutputStreamSettings.
     * @param stderr    The capture of stderr as defined by the stderr OutputStreamSettings.
     */
    public ProcessOutput(int exitValue, StreamOutput stdout, StreamOutput stderr) {
        this.exitValue = exitValue;
        this.stdout = stdout;
        this.stderr = stderr;
    }

    public int getExitValue() {
        return exitValue;
    }

    public StreamOutput getStdout() {
        return stdout;
    }

    public StreamOutput getStderr() {
        return stderr;
    }

    /**
     * Construct a summary message for process output.
     *
     * @param includeStdout whether to include stdout and stderr in summary
     * @return summary message string
     */
    public String getStatusSummary(boolean includeStdout) {
        final StringBuilder message = new StringBuilder();

        if (getExitValue() == 137) {
            // process received SIGKILL, which might indicate OOM
            message.append("\nThe exit code indicates that the process was terminated." +
                    " This may mean the process requires additional memory.\n");
        }

        if (includeStdout) {
            message.append(String.format("\nStdout: %s\nStderr: %s",
                    getStdout().getBufferString(),
                    getStderr().getBufferString()));
        }

        return message.toString();
    }
}
