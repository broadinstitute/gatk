package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.table.TableCodec;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;

import java.util.Arrays;

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

    public XsvLocatableTableCodec() {
        AsciiFeatureCodec<XsvTableFeature>(XsvTableFeature.class);
    }

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
    protected TableFeature createTableFeatureFromSplitLine(final String[] splitLine) {
        return new TableFeature(new SimpleInterval(splitLine[0]), Arrays.asList(splitLine), header);
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
