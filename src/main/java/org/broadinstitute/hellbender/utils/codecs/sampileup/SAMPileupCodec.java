package org.broadinstitute.hellbender.utils.codecs.sampileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.*;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Decoder for single sample SAM pileup data.
 *
 * See the <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format documentation</a> for more details
 * on the format
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SAMPileupCodec extends AsciiFeatureCodec<SAMPileupFeature> {

    // including also space as delimiter
    private static final Pattern SPLIT_PATTERN = Pattern.compile("\\t|( +)");

    // although the minimum number of fields is 6 (see specifications), if the coverage is 0 sometimes bases and qualities does not appear
    private static int MINIMUM_FIELDS = 4;
    // number of maximum fields expected
    private static final int MAXIMUM_FIELDS = 6;

    // regular expression for indels
    private static final Pattern INDEL_REGEXP = Pattern.compile("([0-9]+).*");
    // codec file extensions
    public static final List<String> SAM_PILEUP_FILE_EXTENSIONS = Arrays.asList("pileup", "mpileup");

    public SAMPileupCodec() {
        super(SAMPileupFeature.class);
    }

    @Override
    public Object readActualHeader(final LineIterator lineIterator) {
        // No header for this format
        return null;
    }

    /**
     * Only files with {@link #SAM_PILEUP_FILE_EXTENSIONS} could be parsed
     * @param path path the file to test for parsability with this codec
     * @return {@code true} if the path has the correct file extension {@link #SAM_PILEUP_FILE_EXTENSIONS}, {@code false} otherwise
     */
    @Override
    public boolean canDecode(final String path) {
        final String noBlockCompressedPath;
        if (IOUtil.hasBlockCompressedExtension(path)) {
            noBlockCompressedPath = FilenameUtils.removeExtension(path).toLowerCase();
        } else {
            noBlockCompressedPath = path.toLowerCase();
        }
        return SAM_PILEUP_FILE_EXTENSIONS.stream().anyMatch(ext -> noBlockCompressedPath.endsWith("."+ext));
    }

    @Override
    public SAMPileupFeature decode(String line) {
        // Split the line
        final String[] tokens = SPLIT_PATTERN.split(line.trim(), -1);
        // check the number of fields
        if(tokens.length < MINIMUM_FIELDS || tokens.length > MAXIMUM_FIELDS) {
            throw new CodecLineParsingException(String.format("The SAM pileup line didn't have the expected number of columns (%s-%s): %s. Note that this codes is only valid for single-sample pileups",
                    MINIMUM_FIELDS, MAXIMUM_FIELDS, line));
        }
        // starting parsing
        final String chr = tokens[0];
        final int pos = parseInteger(tokens[1], "position");
        final byte ref = parseBase(tokens[2], "reference");
        final int cov = parseInteger(tokens[3], "coverage");
        // we end parsing here if coverage is 0
        if(cov == 0) {
            return new SAMPileupFeature(chr, pos, ref, new ArrayList<>());
        }
        // parse the elements
        final List<SAMPileupElement> pileupElements = parseBasesAndQuals(tokens[4], tokens[5], ref);
        if(cov != pileupElements.size()) {
            throw new CodecLineParsingException("THe SAM pileup line didn't have the same number of elements as the expected coverage: " + cov);
        }
        return new SAMPileupFeature(tokens[0], pos, ref, pileupElements);

    }

    @VisibleForTesting
    protected List<SAMPileupElement> parseBasesAndQuals(final String bases, final String qualities, final byte ref) {
        try {
            final List<SAMPileupElement> pileupElements = new ArrayList<>(qualities.length());
            int i = 0, j = 0;
            for (; i < bases.length(); i++) {
                final char c = bases.charAt(i);
                switch (c) {
                    case '$': // end of read
                        break;
                    case '^': // mapping quality
                        i++; // the next symbol is a quality value
                        break;
                    case '.':
                    case ',':   // matches reference
                        pileupElements.add(new SAMPileupElement(ref, (byte) SAMUtils.fastqToPhred(qualities.charAt(j++))));
                        break;
                    case '*': // deletion placeholder
                        pileupElements.add(new SAMPileupElement(BaseUtils.Base.D.base, (byte) SAMUtils.fastqToPhred(qualities.charAt(j++))));
                        break;
                    case '+':
                    case '-': // start of indel
                        final String rest = bases.substring(i + 1);
                        final Matcher match = INDEL_REGEXP.matcher(rest);
                        if (!match.matches()) {
                            throw new CodecLineParsingException("The SAM pileup line has an indel marker (+/-) without length");
                        }
                        final String lengthString = match.group(1);
                        final int indelLength = parseInteger(lengthString, "indel-length");
                        i += indelLength + lengthString.length(); // length of number + that many bases + +/- at the start (included in the next i++)
                        break;
                    default:
                        final byte base = parseBase((byte) bases.charAt(i), "reads String");
                        pileupElements.add(new SAMPileupElement(base, (byte) SAMUtils.fastqToPhred(qualities.charAt(j++))));
                }
            }
            if (i != bases.length() || j != qualities.length()) {
                throw new CodecLineParsingException("Not all bases/qualities have been parsed because of a malformed line");
            }
            return pileupElements;
        } catch(IndexOutOfBoundsException e) {
            throw new CodecLineParsingException("Malformed SAM pileup: Different number of bases and qualities found.");
        }
    }

    /**
     * For fast indexing
     */
    @Override
    public Feature decodeLoc(final LineIterator lineIterator) throws IOException {
        String[] tokens = SPLIT_PATTERN.split(lineIterator.next(), -1);
        final int pos = parseInteger(tokens[1], "position");
        return new SimpleFeature(tokens[0], pos, pos);

    }

    private byte parseBase(final byte base, final String parsedValue) {
        if(BaseUtils.isNBase(base)) {
            return BaseUtils.Base.N.base;
        }
        final int index = BaseUtils.simpleBaseToBaseIndex(base);
        if(index == -1) {
            throw new CodecLineParsingException("The SAM pileup line had wrong base at " + parsedValue + ": " + (char) base);
        }
        return BaseUtils.baseIndexToSimpleBase(index);
    }

    private int parseInteger(final String token, final String parsedValue) {
        try {
            return Integer.parseInt(token);
        } catch(NumberFormatException e) {
            throw new CodecLineParsingException("The SAM pileup line had unexpected " + parsedValue + ": " + token);
        }
    }

    private byte parseBase(final String token, final String parsedValue) {
        if(token.length() != 1)  {
            throw new CodecLineParsingException("The SAM pileup line had unexpected base at " + parsedValue + ": " + token);
        }
        return parseBase((byte) token.charAt(0), parsedValue);
    }

    @Override
    public TabixFormat getTabixFormat() {
        return new TabixFormat(TabixFormat.GENERIC_FLAGS, 1, 2, 0, '#', 0);
    }

}
