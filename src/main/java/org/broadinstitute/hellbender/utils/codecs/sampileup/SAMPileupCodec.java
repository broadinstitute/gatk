package org.broadinstitute.hellbender.utils.codecs.sampileup;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.ParsingUtils;

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
 * @author Matt Hanna, Geraldine VdAuwera
 * @since 2014
 */
public class SAMPileupCodec extends AsciiFeatureCodec<SAMPileupFeature> {
    // number of tokens expected
    private static final int TOKEN_COUNTS = 6;
    // delimiter for pileup
    private static final char DELIMITER = '\t';
    // regular expression for indels
    private static final Pattern INDEL_REGEXP = Pattern.compile("([0-9]+).*");
    // codec file extensions
    // TODO: remove .samp
    protected static final List<String> FILE_EXTENSIONS = Arrays.asList("pileup", "mpileup", "samp");

    public SAMPileupCodec() {
        super(SAMPileupFeature.class);
    }

    public SAMPileupFeature decode(String line) {
        //+1 because we want to know if we have more than the max
        final String[] tokens = new String[TOKEN_COUNTS+1];
        // split the line
        final int count = ParsingUtils.split(line, tokens, DELIMITER);
        if(count != TOKEN_COUNTS) {
            throw new CodecLineParsingException("The SAM pileup line didn't have the expected number of tokens. " +
                    "Expected = " + TOKEN_COUNTS + " for basic pileup, but found " + count + ". Line: " + line);
        }
        // create the pileup feature
        final SAMPileupFeature feature = new SAMPileupFeature();
        // set the chromosome and the position
        feature.setChr(tokens[0]);
        // TODO: catch exception parsing the integer
        feature.setStart(Integer.parseInt(tokens[1]));
        // check the reference base
        if(tokens[2].length() != 1)  {
            // TODO: check also the character of the reference to be sure that it is valid
            throw new CodecLineParsingException("The SAM pileup line had unexpected base " + tokens[2] + " on line = " + line);
        }
        feature.setRef(tokens[2].charAt(0));
        parseBasesAndQuals(feature, tokens[4], tokens[5]);
        return feature;
    }

    /**
     * Can the file be decoded?
     * @param path path the file to test for parsability with this codec
     * @return true if the path has the correct file extension, false otherwise
     */
    @Override
    public boolean canDecode(final String path) {
        return FILE_EXTENSIONS.stream().anyMatch(ext -> path.endsWith("."+ext));
    }

    @Override
    public Object readActualHeader(final LineIterator lineIterator) {
        // No header for this format
        return null;
    }

    // TODO: parse base and qualities as PileupElements
    private void parseBasesAndQuals(final SAMPileupFeature feature, final String bases, final String quals) {
        // needs to convert the base string with its . and , to the ref base
        final StringBuilder baseBuilder = new StringBuilder();
        final StringBuilder qualBuilder = new StringBuilder();
        boolean done = false;
        for ( int i = 0, j = 0; i < bases.length() && ! done; i++ ) {
            final char c = bases.charAt(i);
            switch ( c ) {
                case '.':   // matches reference
                case ',':   // matches reference
                    baseBuilder.append(feature.getRef());
                    qualBuilder.append(quals.charAt(j++));
                    break;
                case '$':   // end of read
                    break;
                case '*':   // end of indel?
                    j++;
                    break;
                case '^':   // mapping quality
                    i++;
                    break;
                case '+':   // start of indel
                case '-':   // start of indel
                    // match the expression for indels
                    final String rest = bases.substring(i+1);
                    Matcher match = INDEL_REGEXP.matcher(rest);
                    if ( ! match.matches() ) {
                        // TODO: is this correct?
                        if ( feature.getRef() != '*' ) {
                            throw new CodecLineParsingException("Bad pileup format: " + bases + " at position " + i);
                        }
                        // TODO: is this really done? if so, why it is needed the flag and not a break?
                        done = true;
                    } else {
                        String g = match.group(1);
                        int l = Integer.parseInt(g);
                        i += l + g.length();    // length of number + that many bases + +/- at the start (included in the next i++)
                    }
                    break;
                default:   // non reference base
                    baseBuilder.append(c);
                    qualBuilder.append(quals.charAt(j++));
            }
        }
        feature.setPileupBases(baseBuilder.toString());
        feature.setPileupQuals(qualBuilder.toString());
    }
}
