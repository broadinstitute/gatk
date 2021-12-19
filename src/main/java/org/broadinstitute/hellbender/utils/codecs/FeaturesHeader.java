package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

public class FeaturesHeader {
    private final String className;
    private final String version;
    private final SAMSequenceDictionary dictionary;
    private final List<String> sampleNames;

    public FeaturesHeader( final String className,
                           final String version,
                           final SAMSequenceDictionary dictionary,
                           final List<String> sampleNames ) {
        Utils.nonNull(className);
        Utils.nonNull(version);
        this.className = className;
        this.version = version;
        this.dictionary = dictionary;
        this.sampleNames = sampleNames;
    }

    public String getClassName() {
        return className;
    }

    public String getVersion() {
        return version;
    }

    public SAMSequenceDictionary getDictionary() {
        return dictionary;
    }

    public List<String> getSampleNames() {
        return sampleNames;
    }
}
