package org.broadinstitute.hellbender.utils.codecs.gtf;

import java.util.LinkedHashMap;
import java.util.Map;

public class Gff3Codec extends GffCodec {
    private final static String GFF3_FILE_EXTENSION = "gff";

    public Gff3Codec() {
        super(GFF3_FILE_EXTENSION);
    }

    @Override
    public Map<String,String> parseAttributes(final String attributesString) {
        final Map<String, String> attributes = new LinkedHashMap<>();
        final String[] splitLine = attributesString.split(";");
        for(String attribute : splitLine) {
            final String[] key_value = attribute.split("=");
            if (key_value.length<2) {
                continue;
            }
            attributes.put(key_value[0], key_value[1]);
        }
        return attributes;
    }

    @Override
    public String getFirstLineStart() {
        return "##gff-version 3";
    }
}
