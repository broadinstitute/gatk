package org.broadinstitute.hellbender.fakedata;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class SimulatedSamples {
    public static List<String> phonySamples(final int numSamples) {
        return IntStream.range(0, numSamples)
                .mapToObj(n -> String.format("Sample%d", n)).collect(Collectors.toList());
    }
}
