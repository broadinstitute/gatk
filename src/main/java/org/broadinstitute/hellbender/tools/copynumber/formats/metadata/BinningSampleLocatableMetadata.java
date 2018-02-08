package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountType;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;

import java.util.List;

/**
 * Interface for marking objects that contain metadata associated with a collection of locatables
 * and with a single sample, as well as read count type and list of binning configurations
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public interface BinningSampleLocatableMetadata extends SampleLocatableMetadata {

    /**
     * @return list of binning configurations
     */
    List<ReadCountCovariateBinningConfiguration> getCovariateBinningConfigurations();

    /**
     * @return read count type
     */
    ReadCountType getReadCountType();
}
