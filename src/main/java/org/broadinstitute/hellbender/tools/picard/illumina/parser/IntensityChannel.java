package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * The channels in a FourChannelIntensityData object, and the channels produced by a ClusterIntensityFileReader,
 * for cases in which it is desirable to handle these abstractly rather than having the specific names
 * in the source code.
 *
 * @author alecw@broadinstitute.org
 */
public enum IntensityChannel {

    A, C, G, T;

    public static final int NUM_CHANNELS = values().length;
}
