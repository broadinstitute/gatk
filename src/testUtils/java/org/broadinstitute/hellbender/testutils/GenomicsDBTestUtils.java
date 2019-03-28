package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.Collections;
import java.util.List;


public final class GenomicsDBTestUtils {

    /**
     * don't instantiate a utility class
     */
    private GenomicsDBTestUtils(){}

    /**
     * @param workspace the workspace folder of a GenomicsDB with a single array and standard json file layout
     * @return a string formatted as a genomicsDB uri pointing to the given workspace i.e "gendb:///pathTo/workspace
     */
    public static String makeGenomicsDBUri(final File workspace){
        return IOUtils.GENOMIC_DB_URI_SCHEME + "://" + workspace.getAbsolutePath();
    }

    /**
     * Create a temporary GenomicsDB containing a single interval of data from a set of gvcfs
     * this database will be deleted on jvm shutdown automatically
     * @param gvcfs, a gvcf to load from
     * @param interval the interval to load
     * @return the created workspace folder containing the new GenomicsDB
     */
    public static File createTempGenomicsDB(final File gvcfs, final Locatable interval) {
        return createTempGenomicsDB(Collections.singletonList(gvcfs), interval);
    }

    /**
     * Create a temporary GenomicsDB containing a single interval of data from a set of gvcfs
     * this database will be deleted on jvm shutdown automatically
     * @param gvcfs, a List of a GVCFs to load from
     * @param interval the interval to load
     * @return the created workspace folder containing the new GenomicsDB
     */
    public static File createTempGenomicsDB(final List<File> gvcfs, final Locatable interval) {
        final File workspaceDir = BaseTest.createTempDir("genomicsDBWorkspace");

        final CommandLineProgramTester importer = GenomicsDBImport.class::getSimpleName;

        final ArgumentsBuilder args = new ArgumentsBuilder();
        gvcfs.forEach(args::addVCF);


        final String workspace = new File(workspaceDir, "workspace").getAbsolutePath();
        args.addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);
        args.addInterval(interval);
        importer.runCommandLine(args);
        return new File(workspace);
    }

    /**
     * Create a temporary GenomicsDB containing multiple intervals of data from a set of gvcfs
     * this database will be deleted on jvm shutdown automatically
     * @param gvcf, a GVCF to load from
     * @param intervals the list of intervals to load
     * @param mergeIntervals if true, import data over the span of all intervals
     * @return the created workspace folder containing the new GenomicsDB
     */
    public static File createTempGenomicsDB(final File gvcf, final List<Locatable> intervals, final boolean mergeIntervals) {
        final File workspaceDir = BaseTest.createTempDir("genomicsDBWorkspace");

        final CommandLineProgramTester importer = GenomicsDBImport.class::getSimpleName;

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addVCF(gvcf);

        final String workspace = new File(workspaceDir, "workspace").getAbsolutePath();
        args.addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);
        intervals.forEach(args::addInterval);
        if (mergeIntervals) {
            args.addBooleanArgument(GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME, true);
        }
        importer.runCommandLine(args);
        return new File(workspace);
    }
}
