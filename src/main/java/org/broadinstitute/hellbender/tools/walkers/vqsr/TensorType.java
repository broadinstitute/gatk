package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.barclay.argparser.CommandLineParser;

/**
 * TensorType documents the tensors available and what information they encode.
 * Created by sam on 2/26/18.
 */
public enum TensorType implements CommandLineParser.ClpEnum {
    reference(4, false, " 1 Hot encoding of a reference sequence. "),
    read_tensor(15, true, "Read tensor are 3D tensors spanning aligned reads, sites and channels. " +
            "The maximum number of reads is a hyper-parameter typically set to 128. " +
            "There are 15 channels in the read tensor. They correspond to the reference sequence data (4)," +
            " read sequence data (4), insertions and deletions (2) read flags (4) and mapping quality (1).")
    ;

    private int channels;
    private boolean readsRequired;
    private String description;

    /** The tensor type documents how and what information is encoded in numeric tensors for processing
     * by neural networks.
     *
     * @param channels The number of channels in this tensor type. For example,
     *                 1-hot encoded DNA will have 4 channels, one for each nucleotide base.
     *                 A standard image has 3 channels for red, green and blue.
     * @param readsRequired True if this tensor encodes read-level data from a SAM/BAM, false otherwise.
     * @param description Documentation describing the tensor, should explain what each of the channels mean.
     */
    TensorType(int channels, boolean readsRequired, String description) {
        this.channels = channels;
        this.readsRequired = readsRequired;
        this.description = description;
    }

    @Override
    public String getHelpDoc() {
        return description;
    }

    public boolean isReadsRequired() {
        return readsRequired;
    }

    public int getChannels() { return channels; }

}
