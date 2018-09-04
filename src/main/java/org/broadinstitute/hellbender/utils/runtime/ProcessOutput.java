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
}
