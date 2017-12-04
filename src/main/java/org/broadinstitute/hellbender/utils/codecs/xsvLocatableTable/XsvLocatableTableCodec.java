package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;

/**
 * Created by jonn on 12/4/17.
 */
public final class XsvLocatableTableCodec extends TableCodec<XsvTableFeature> {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public boolean canDecode(final String path) {
        // Check for a sibling config file with the same name, .config as extension
        // Open that config file
        // Validate config file
        //     Expected keys present
        //     Key values are valid
        return false;
    }

    @Override
    public XsvTableFeature decode(final String s) {
        return null;
    }

    @Override
    public Object readActualHeader(final LineIterator reader) {
        return null;
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
