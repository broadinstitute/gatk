package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.CommandLineParser;

/**
 * Created by sam on 2/26/18.
 */
public enum TensorMapEnum implements CommandLineParser.ClpEnum {
    reference(4, " 1 Hot encoding of a reference sequence. "),
    read_tensor(15, "Read tensor are 3D tensors spanning aligned reads, sites and channels. " +
            "The maximum number of reads is a hyper- parameter typically set to 128. " +
            "There are 17 channels in the read tensor. They correspond to the reference sequence data (4)," +
            " read sequence data (4), insertions and deletions (2) read flags (4) and mapping quality (1).")
    ;

    private int channels;
    private String description;
    TensorMapEnum(int channels, String description) {
        this.channels = channels;
        this.description = description;
    }



    @Override
    public String getHelpDoc() {
        return description;
    }


}
