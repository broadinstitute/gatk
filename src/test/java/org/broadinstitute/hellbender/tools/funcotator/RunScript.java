//package org.broadinstitute.hellbender.tools.funcotator;
//
////import org.apache.commons.exec.CommanLine;
////import org.apache.commons.exec.DefaultExecutor;
////import org.apache.commons.exec.ExecuteException;
//import java.io.BufferedReader;
//import java.io.IOException;
//import java.io.InputStreamReader;
//import sun.tools.jar.CommandLine;
//import java.io.IOException;
//import org.omg.CORBA.portable.InputStream;
//
//public class RunScript {
//    int exitValue;
//    String commandString;
//
//    public void runScript(String command) {
//        commandString = command;
//        CommandLine cmdLine = CommandLine.parse(commandString);
//        DefaultExecutor defaultExecutor = new DefaultExecutor();
//        defaultExecutor.setExitValue(0);
//        try {
//            exitValue = defaultExecutor.execute(cmdLine);
//        } catch (ExecuteException e) {
//            System.err.println("Execution failed.");
//            e.printStackTrace();
//        } catch (IOException e) {
//            System.err.println("permission denied.");
//            e.printStackTrace();
//        }
//    }
//}
