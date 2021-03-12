package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SamFiles;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadsBundle;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;

/**
 * Create a JSON bundle file for use with GATK tools.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Create a reads JSON bundle files for use with GATK tool",
        oneLineSummary = "Create reads JSON bundle files for use with GATK tools",
        programGroup = OtherProgramGroup.class
)
public class CreateReadsBundleJson extends CommandLineProgram {
    public static final String NO_INDEX_FULL_NAME = "no-index";

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="Path to SAM/BAM/CRAM to create a reads-bundle.json for")
    GATKPath reads;

    @Argument(fullName = StandardArgumentDefinitions.READ_INDEX_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_INDEX_SHORT_NAME,
            doc = "Path to index of BAM/CRAM specified with " + StandardArgumentDefinitions.INPUT_LONG_NAME
                    + ". If not specified the index will be automatically inferred.",
            optional = true,
            mutex = {"no-index"})
    GATKPath index;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Path the output bundle file.")
    GATKPath outputBundle;

    @Argument(fullName = "no-index",
            doc =" this must be specified to create a bundle without an index",
            optional = true,
            mutex = {StandardArgumentDefinitions.READ_INDEX_LONG_NAME})
    boolean noIndex = false;

    @Override
    protected Object doWork() {
        if( index == null && !noIndex){
            index = IOUtils.toGATKPath(SamFiles.findIndex(reads.toPath()));
            if (index == null){
                throw new UserException.MissingIndex("Could not locate an index for " + reads.getRawInputString() +
                        ", either specify the index path with --" +StandardArgumentDefinitions.READ_INDEX_LONG_NAME + " or specify the "
                        + "--" + NO_INDEX_FULL_NAME + " argument.");
            }
        }

        //TODO: should the UT8-encoding be a constant that lives somewhere else ?
        try (final OutputStream os = outputBundle.getOutputStream();
             final OutputStreamWriter sw = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
            final ReadsBundle bundle = new ReadsBundle(reads, index);
            sw.write(bundle.toJSON());
        } catch (final IOException e) {
            throw new UserException(String.format("Failed writing bundle to output %s", outputBundle), e);
        }

        return null;
    }
}
