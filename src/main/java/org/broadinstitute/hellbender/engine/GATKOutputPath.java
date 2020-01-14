package org.broadinstitute.hellbender.engine;

/**
 * Marker class used to distinguish command line arguments that are tool inputs from tool outputs.
 */
public class GATKOutputPath extends GATKPathSpecifier {
    private static final long serialVersionUID = 1L;

    /**
     * @param uriString The raw value for this output path as provided by the user. Can be a local file or valid URI.
     */
    public GATKOutputPath(final String uriString) {
        super(uriString);
    }
    
}
