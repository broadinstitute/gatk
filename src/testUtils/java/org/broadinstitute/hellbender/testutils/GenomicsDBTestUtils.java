package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.IntervalUtils;

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
        return FeatureDataSource.GENOMIC_DB_URI_SCHEME + workspace.getAbsolutePath();
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
        args.addArgument("L", IntervalUtils.locatableToString(interval));
        importer.runCommandLine(args);
        return new File(workspace);
    }
}
