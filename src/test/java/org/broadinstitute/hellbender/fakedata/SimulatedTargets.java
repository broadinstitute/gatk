package org.broadinstitute.hellbender.fakedata;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class SimulatedTargets {

    // phony targets whose intervals don't really matter
    public static List<Target> phonyTargets(final int numTargets) {
        return IntStream.range(0, numTargets).mapToObj(n -> {
            final String name = String.format("Target%d", n);
            final SimpleInterval interval = new SimpleInterval("chr", 2*n + 1, 2*n + 2);
            return new Target(name, interval);
        }).collect(Collectors.toList());
    }
}
