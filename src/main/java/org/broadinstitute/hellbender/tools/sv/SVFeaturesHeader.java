package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class SVFeaturesHeader {
    private final String className;
    private final String version;
    private final SAMSequenceDictionary dictionary;
    private final List<String> sampleNames;
    private Map<String, Integer> sampleNameMap;

    public SVFeaturesHeader( final String className,
                             final String version,
                             final SAMSequenceDictionary dictionary,
                             final List<String> sampleNames ) {
        Utils.nonNull(className);
        Utils.nonNull(version);
        this.className = className;
        this.version = version;
        this.dictionary = dictionary;
        this.sampleNames = new ArrayList<>(sampleNames);
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

    public Integer getSampleIndex( final String sampleName ) {
        if ( sampleNameMap == null ) {
            int nNames = sampleNames.size();
            sampleNameMap = new HashMap<>((nNames * 4 + 2) / 3);
            for ( int idx = 0; idx != nNames; ++idx ) {
                sampleNameMap.put(sampleNames.get(idx), idx);
            }
        }
        return sampleNameMap.get(sampleName);
    }
}
