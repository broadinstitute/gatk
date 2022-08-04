package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MummerExecutor {

    private File mummerExecutableDirectory;

    public MummerExecutor(File mummerExecutableDirectory){
        this.mummerExecutableDirectory = mummerExecutableDirectory;
    }

    public Path executeMummer(File fasta1, File fasta2){
        // runShellCommand(dnadiff) -- might need to add perl if doesn't work
        return null;
    }

    public static void runShellCommand(String[] command, Map<String, String> environment, File stdoutCaptureFile, boolean printStdout){
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
    }

    public static void runPythonCommand(String script, List<String> scriptArguments, Map<String, String> additionalEnvironmentVars, File stdoutCaptureFile, boolean printStdout){
        Map<String, String> environment = new HashMap<>();
        environment.putAll(System.getenv());
        if(additionalEnvironmentVars != null){
            environment.putAll(additionalEnvironmentVars);
        }
        List<String> args = new ArrayList<>();
        args.add("python");
        args.add(script);
        args.addAll(scriptArguments);
        runShellCommand(args.toArray(new String[]{}), environment, stdoutCaptureFile, printStdout);
        /*PythonScriptExecutor executor = new PythonScriptExecutor(true);
        boolean status = executor.executeScript(script, null, Arrays.asList(scriptArguments));*/
    }

}
