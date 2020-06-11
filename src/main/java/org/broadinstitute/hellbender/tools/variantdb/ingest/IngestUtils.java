package org.broadinstitute.hellbender.tools.variantdb.ingest;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;

public class IngestUtils {
    static final Logger logger = LogManager.getLogger(IngestUtils.class);

    public static String getSampleName(final VCFHeader inputVCFHeader) {
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        if (samples.numberOfSamples() > 1) {
            throw new UserException("This tool can only be run on single sample vcfs");
        }
        return samples.getSample(0);
    }

    public static String getSampleId(final String sampleName, final File sampleMap) {
        String sampleId = null;
        //  Because BigQuery only supports partitioning based on timestamp or integer,
        // sample names will be remapped into sample_id integers
        try {
            BufferedReader br = new BufferedReader(new FileReader(sampleMap));

            String line; // Reading header, Ignoring
            while ((line = br.readLine()) != null && !line.isEmpty()) {
                String[] fields = line.split(",");
                String name = fields[1];
                if (sampleName.equals(name)) {
                    sampleId = fields[0];
                    break;
                }
            }
            br.close();
            if (sampleId == null) {
                // sampleName not found
                throw new UserException("Sample " + sampleName + " could not be found in sample mapping file");
            }
        } catch (final IOException ioe) { // FileNotFoundException e,
            throw new UserException("Could not find sample mapping file");
        }
        return sampleId;
    }

    public static Path createSampleDirectory(Path parentDirectory, int sampleDirectoryNumber) {
        // If this sample set directory doesn't exist yet -- create it
        final String sampleDirectoryName = String.valueOf(sampleDirectoryNumber);
        final Path sampleDirectoryPath = parentDirectory.resolve(sampleDirectoryName);
        final File sampleDirectory = new File(sampleDirectoryPath.toString());
        if (!sampleDirectory.exists()) {
            sampleDirectory.mkdir();
        }
        return sampleDirectoryPath;
    }

    // To determine which directory (and ultimately table) the sample's data will go into
    // Since tables have a limited number of samples (default is 4k)
    public static int getTableNumber(String sampleId, int sampleMod) { // this is based on sample id
        // sample ids 1-4000 will go in directory 001
        int sampleIdInt = Integer.valueOf(sampleId); // TODO--should sampleId just get refactored as a long?
        // subtract 1 from the sample id to make it 1-index (or do we want to 0-index?) and add 1 to the dir
        int directoryNumber = Math.floorDiv((sampleIdInt - 1), sampleMod) + 1; // TODO omg write some unit tests
        return directoryNumber;
    }


}
