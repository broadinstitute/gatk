package org.broadinstitute.hellbender.tools.spark.sv.sga;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents an SGA module that can be called via "run" (in base) to do actual work.
 */
final class SGAModule extends ExternalCommandlineProgramModule {

    private final String moduleName;

    public SGAModule(final String moduleName){
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
