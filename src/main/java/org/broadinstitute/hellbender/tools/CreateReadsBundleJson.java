package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SamFiles;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.bundle.ReadsBundle;
import org.broadinstitute.hellbender.utils.io.IOUtils;

public class CreateReadsBundleJson extends CommandLineProgram {
    public static final String NO_INDEX_FULL_NAME = "no-index";

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="Path to SAM/BAM/CRAM to create a reads-bundle.json for",
            optional = false)
    GATKPath reads;

    @Argument(fullName = StandardArgumentDefinitions.READ_INDEX_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_INDEX_SHORT_NAME,
            doc = "Path to index of BAM/CRAM specified with " + StandardArgumentDefinitions.INPUT_LONG_NAME
                    + ". If not specified the index will be automatically inferred.",
    mutex = {"no-index"})
    GATKPath index;

    @Argument(fullName = "no-index", doc =" this must be specified to create a bundle without an index", mutex = {StandardArgumentDefinitions.READ_INDEX_LONG_NAME})
    boolean noIndex = false;

    @Override
    protected Object doWork() {
        final ReadsBundle bundle;
        if( index == null && !noIndex){
            index = IOUtils.toGATKPath(SamFiles.findIndex(reads.toPath()));
            if (index == null){
                throw new UserException.MissingIndex("Could not locate an index for " + reads.getRawInputString() +
                        ", either specify the index path with --" +StandardArgumentDefinitions.READ_INDEX_LONG_NAME + " or specify the "
                        + "--" + NO_INDEX_FULL_NAME + " argument.");
            }
        }
        return null;
    }
}
