package org.broadinstitute.hellbender.utils.codecs.gtf;

import java.util.HashMap;
import java.util.Map;

public class Gff3Codec extends GffCodec {

    @Override
    public Map<String,String> parseAttributes(final String attributesString) {
        final Map<String, String> attributes = new HashMap<>();
        final String[] splitLine = attributesString.split(";");
        for(String attribute : splitLine) {
            final String[] key_value = attribute.split("=");
            attributes.put(key_value[0], key_value[1]);
        }
        return attributes;
    }

    @Override
    public String getFirstLineStart() {
        return "#gff-version 3";
    }
}
