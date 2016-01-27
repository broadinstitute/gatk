package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

/**
 * Read XHMM discovery output formatted file.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class XHMMOutputReader extends TableReader<XHMMOutputRecord> {

    public XHMMOutputReader(File file) throws IOException {
        super(file);
    }

    @Override
    protected XHMMOutputRecord createRecord(final DataLine dataLine) {
        return new XHMMOutputRecord(dataLine);
    }
}
