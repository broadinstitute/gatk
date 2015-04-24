package org.broadinstitute.hellbender.utils.report;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.commandline.Gatherer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

public final class GATKReportGatherer extends Gatherer {
    @Override
    public void gather(List<File> inputs, File output) {
        //Combines inputs GATKReport to one output

        PrintStream o;
        try {
            o = new PrintStream(output);
        } catch (FileNotFoundException e) {
            throw new UserException(String.format("File %s to be output by GATKReportGatherer function was not found", output));
        }

        GATKReport current = new GATKReport();
        boolean isFirst = true;
        for (File input : inputs) {
            if (isFirst) {
                current = new GATKReport(input);
                isFirst = false;
            } else {
                current.concat(new GATKReport(input));
            }
        }

        current.print(o);
        o.close();
    }
}
