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

public class IngestUtils {
    static final Logger logger = LogManager.getLogger(IngestUtils.class);

    public static String getSampleName(final VCFHeader inputVCFHeader) {
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        if (samples.numberOfSamples() > 1){
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


}
