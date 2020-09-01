package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceDownloader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

/**
 * Class to hold integration tests for the dbSNP data source.
 */
public class DbSnpIntegrationTest extends CommandLineProgramTest {

    private final int BUFFER_SIZE = 1024 * 1024;

    private final String DB_SNP_FILE_NAME = "hg38_All_20180418.vcf.gz";
    private final String DB_SNP_INDEX_FILE_NAME = "hg38_All_20180418.vcf.gz.tbi";

    /**
     * Extract the dbSNP VCF file (and the associated index) from the funcotator datasources and return the path to the
     * extracted file.
     * Adapted from a stack-overflow post:
     * https://stackoverflow.com/a/41853978
     *
     * Note: We return the path to the dbSNP file, but the index will also be extracted in parallel so that GATK can
     * properly parse the file.
     *
     * @param in {@link InputStream} pointing to Funcotator Data Sources on gcloud.
     * @return A {@link Path} object pointing to the dbSNP VCF extracted from the given {@link InputStream}.
     */
    private Path extractAndReturnDbSnpPath(final InputStream in) {

        boolean foundDbSnp = false;
        Path dbSnpPath = null;

        logger.info("Getting info from gcloud now...");

        try {
            final GzipCompressorInputStream gzipInputStream = new GzipCompressorInputStream(in);
            try ( final TarArchiveInputStream tarInputStream = new TarArchiveInputStream(gzipInputStream) ) {
                TarArchiveEntry entry = (TarArchiveEntry) tarInputStream.getNextEntry();

                logger.info("Connected to gcloud archive file.");

                while ( entry != null ) {

                    logger.info("Current archive file: " + entry.getName());

//                    if ( entry.isDirectory() ) {
//                        final File    f       = new File(entry.getName());
//                        final boolean created = f.mkdir();
//                        if ( !created ) {
//                            System.out.printf("Unable to create directory '%s', during extraction of archive contents.\n",
//                                    f.getAbsolutePath());
//                        }
//                    }
//                    else {
//                        int count;
//
//                        final byte[] data = new byte[ BUFFER_SIZE ];
//                        final FileOutputStream fos  = new FileOutputStream(entry.getName(), false);
//                        try ( final BufferedOutputStream dest = new BufferedOutputStream(fos, BUFFER_SIZE) ) {
//                            while ( (count = tarIn.read(data, 0, BUFFER_SIZE)) != -1 ) {
//                                dest.write(data, 0, count);
//                            }
//                        }
//                    }
                    // Get the next file in the archive:
                    entry = (TarArchiveEntry) tarInputStream.getNextEntry();
                }
            }
        }
        catch (final IOException ex) {
            throw new UserException("Unable to read contents of tgz file.", ex);
        }

        if ( !foundDbSnp ) {
            throw new UserException("Unable to find dbSNP data source in the data sources tar.gz file!");
        }

        return dbSnpPath;
    }

    @Test(groups={"cloud"})
//    @Test
    public void testDbSnpDataSourceParsing() {
        // 1 - Extract the dbSNP file from the current datasources for Funcotator:
        logger.info("Creating input stream from gcloud file:");
        try (final InputStream dataSourcesInputStream = new BufferedInputStream(Files.newInputStream(FuncotatorDataSourceDownloader.SOMATIC_GCLOUD_DATASOURCES_PATH))) {
            final Path dbSnpPath = extractAndReturnDbSnpPath(dataSourcesInputStream);

            // 2 - Create a FeatureDataSource from the dbSNP VCF:
            final FeatureDataSource<VariantContext> dbSnpDataSource = new FeatureDataSource<>(dbSnpPath.toUri().toString());

            // Do a dummy check here:
            Assert.assertNotNull(dbSnpDataSource);

            // 3 - Attempt to read sites and features from the FeatureDataSource that would fail with the old code:
            final List<Locatable> intervalsToQuery = Arrays.asList(
                    new SimpleInterval("chr1", 84349784, 84349786),
                    new SimpleInterval("chrX", 80688069, 80688071), // rs3122407 T/C
                    new SimpleInterval("chrY", 13355943, 13355945), // rs2032654 A/G
                    new SimpleInterval("chrM", 5131, 5134)    // rs199476116 2-BP DEL, 5132AA
            );

            for (int i = 0; i < intervalsToQuery.size(); ++i) {
                final Locatable interval = intervalsToQuery.get(i);
                final List<VariantContext> features = dbSnpDataSource.queryAndPrefetch(interval);

                Assert.assertEquals(features.size(), 1);
            }
        }
        catch (final IOException ex) {
            throw new UserException("Unable to open data sources from gcloud!", ex);
        }
    }
}
