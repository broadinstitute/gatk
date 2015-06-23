package org.broadinstitute.hellbender.engine.dataflow.datasources;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.nio.charset.Charset;

public class FakeReferenceSource {
    public static ReferenceBases bases(SimpleInterval interval) {
        StringBuilder stringBuilder = new StringBuilder();
        int start = interval.getStart();
        int end = interval.getEnd();
        for (int i = start; i <= end; ++i) {
            stringBuilder.append("A");
        }

        return new ReferenceBases(stringBuilder.toString().getBytes(Charset.forName("UTF-8")), interval);
    }

}
