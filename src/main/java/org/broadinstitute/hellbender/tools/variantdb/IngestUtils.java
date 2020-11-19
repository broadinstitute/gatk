package org.broadinstitute.hellbender.tools.variantdb;

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
                // sampleId not found
                throw new UserException("Sample " + sampleId + " could not be found in sample mapping file");
            }
        } catch (final IOException ioe) { // FileNotFoundException e,
            throw new UserException("Could not find sample mapping file");
        }
        return sampleId;
    }

    // To determine which table the sample's data will go into
    // Since tables have a limited number of samples (default is 4k)
    public static int getTableNumber(String sampleId, int sampleMod) { // this is based on sample id
        // sample ids 1-4000 will go in directory 001
        long sampleIdInt = Long.valueOf(sampleId);
        return getTableNumber(sampleIdInt, sampleMod);
    }

    public static int getTableNumber(long sampleId, int sampleMod) { // this is based on sample id
        // sample ids 1-4000 will go in directory 001
        // subtract 1 from the sample id to make it 1-index (or do we want to 0-index?) and add 1 to the dir
        int directoryNumber = new Long(Math.floorDiv((sampleId - 1), sampleMod) + 1).intValue(); // TODO omg write some unit tests
        return directoryNumber;
    }


}
