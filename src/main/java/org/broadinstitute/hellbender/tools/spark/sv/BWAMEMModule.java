package org.broadinstitute.hellbender.tools.spark.sv;


import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents an BWA-MEM module that can be called via "run" (in base) to do actual alignment work.
 */
public final class BWAMEMModule extends ExternalCommandlineProgramModule {

    @Override
    public List<String> initializeCommands(final Path pathToProgram){
        final ArrayList<String> res = new ArrayList<>();
        res.add(pathToProgram.toString());
        res.add("mem");
        return res;
    }
}