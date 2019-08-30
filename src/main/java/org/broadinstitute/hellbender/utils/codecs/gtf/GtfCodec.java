package org.broadinstitute.hellbender.utils.codecs.gtf;

import com.sun.tools.javah.Gen;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.SimpleFeature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfFeature;
import scala.Int;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GtfCodec extends GffCodec {
    private final static String GTF_FILE_EXTENSION = "gtf";

    public GtfCodec() {
        super(GTF_FILE_EXTENSION);
    }

    @Override
    public Map<String,String> parseAttributes(final String attributesString) {
        final Map<String, String> attributes = new HashMap<>();

        int currentSubstringStart=0;
        while (currentSubstringStart < attributesString.length()) {
            final int keyStart = currentSubstringStart;
            final int keyEnd = attributesString.indexOf(" ");
            final String key = attributesString.substring(keyStart, keyEnd).trim();

            final int valueStart = attributesString.indexOf("\"", keyEnd);
            final int valueEnd = attributesString.indexOf("\"", valueStart + 1);
            final String value = attributesString.substring(valueStart+1, valueEnd).trim();
            attributes.put(key, value);

            currentSubstringStart = valueEnd + 1;

        }
        return attributes;
    }

    @Override
    public String getFirstLineStart() {
        return "#gtf-version 2";
    }
}