package org.broadinstitute.hellbender.utils.illumina;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair;

import java.util.*;

import static htsjdk.samtools.ReservedTagConstants.XT;
import static htsjdk.samtools.util.SequenceUtil.reverseComplement;
import static htsjdk.samtools.util.StringUtil.stringToBytes;
import static java.lang.Integer.valueOf;
import static java.lang.Math.min;
import static java.lang.String.format;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.MIN_MATCH_BASES;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.adapterTrimIlluminaPairedReads;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.adapterTrimIlluminaSingleRead;
import static org.broadinstitute.hellbender.utils.illumina.ClippingUtil.findIndexOfClipSequence;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.DUAL_INDEXED;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.FLUIDIGM;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.INDEXED;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.NEXTERA_V2;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.PAIRED_END;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.SINGLE_END;
import static org.broadinstitute.hellbender.utils.illumina.IlluminaUtil.IlluminaAdapterPair.values;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertNull;
import static org.testng.Assert.assertTrue;

/**
 *
 */
public class ClippingUtilTest {

    @Test(dataProvider = "clipTestData")
    public void testBasicClip(final String testName, final String read, final String clip, final int minMatch, final double errRate, final int expected) {
        final byte[] r = (read == null) ? null : stringToBytes(read);
        final byte[] c = (clip == null) ? null : stringToBytes(clip);

        final int result = findIndexOfClipSequence(r, c, minMatch, errRate);
        assertEquals(result, expected, testName);

    }


    @Test(dataProvider = "clipTestData")
    public void testSingleEndSamRecordClip(final String testName, final String read, final String clip, final int minMatch,
                                           final double errRate, final int expected) {
        if (read == null) return; // Silly case

        final SingleEndAdapter adapter = new SingleEndAdapter(testName, clip);

        for (final boolean reverse : new boolean[]{false, true}) {
            final SAMRecord rec = new SAMRecord(null);
            if (reverse) {
                rec.setReadString(reverseComplement(read));
                rec.setReadNegativeStrandFlag(true);
            } else {
                rec.setReadString(read);
            }

            final AdapterPair matchedAdapter = adapterTrimIlluminaSingleRead(rec, minMatch, errRate, adapter);
            if (expected == -1) {
                assertNull(matchedAdapter, testName);
                assertNull(rec.getAttribute(XT), testName);
            } else {
                assertEquals(matchedAdapter, adapter, testName);
                assertEquals(rec.getAttribute(XT), expected + 1, testName);

                // Test that if minMatch is decreased, it still matches
                for (int i = 1; i < minMatch - 3; ++i) {
                    assertEquals(adapterTrimIlluminaSingleRead(rec, minMatch - i, errRate, adapter), adapter, testName);
                }

                // Skip this test for high error rates, because almost anything will match.
                if (errRate < 0.5) {
                    // Test that if minMatch is increased, it now fails to match
                    assertNull(adapterTrimIlluminaSingleRead(rec, minMatch + 1, errRate, adapter), testName);

                    // Test that if errRate is increased, it still matches
                    assertNotNull(adapterTrimIlluminaSingleRead(rec, minMatch, errRate + 0.1, adapter), testName);
                }

                // Very low error threshold should cause match failure
                assertNull(adapterTrimIlluminaSingleRead(rec, minMatch, 0.01, adapter), testName);
            }
        }

    }

    @Test(dataProvider = "clipPairedTestData")
    public void testPairedEndClip(final String testName, final String read1, final String read2, final AdapterPair expected) {

        final SAMRecord rec1 = new SAMRecord(new SAMFileHeader());
        rec1.setReadString(read1);
        rec1.setFirstOfPairFlag(true);
        final SAMRecord rec2 = new SAMRecord(new SAMFileHeader());
        rec2.setReadString(read2);
        rec2.setSecondOfPairFlag(true);

        final AdapterPair result = adapterTrimIlluminaPairedReads(rec1, rec2,
                INDEXED, PAIRED_END);
        if (result != null) {
            assertEquals(result.getName(), expected.getName(), testName);
        } else {
            assertNull(expected, testName);
        }
    }

    @Test
    public void testOneSidedMatchSupersededByTwoSidedMatch() {
        final int readLength = 36;
        final int adapterLength = 30;
        final String read1 = patchAdapterSubsequenceIntoRead(makeBogusReadString(readLength), SINGLE_END.get3PrimeAdapterInReadOrder(), adapterLength);
        final String read2 = patchAdapterSubsequenceIntoRead(makeBogusReadString(readLength), SINGLE_END.get5PrimeAdapterInReadOrder(), adapterLength);
        final SAMRecord rec1 = new SAMRecord(null);
        rec1.setReadString(read1);
        final SAMRecord rec2 = new SAMRecord(null);
        rec2.setReadString(read2);
        AdapterPair result = adapterTrimIlluminaPairedReads(rec1, rec2, adapterLength - 2, 0.1,
                PAIRED_END, SINGLE_END);
        assertEquals(result, SINGLE_END);

        // Confirm that without SINGLE_END, PAIRED_END would one-sided match, if match length is short enough.
        result = adapterTrimIlluminaPairedReads(rec1, rec2, adapterLength / 2, 0.1, PAIRED_END);
        assertEquals(result, PAIRED_END);

        // However, without shorter match length, one-sided match will fail.
        assertNull(adapterTrimIlluminaPairedReads(rec1, rec2, adapterLength - 2, 0.1, PAIRED_END));
    }

    @Test(dataProvider = "testAdapterInAllReadPositionsDataProvider")
    public void testAdapterInAllReadPositions(final int readLength) {
        final int minAdapterLength = 6;
        for (final IlluminaAdapterPair adapterPair : values()) {
            final AdapterMarker marker = new AdapterMarker(adapterPair);
            for (int adapterPosition = 0; adapterPosition < readLength; ++adapterPosition) {
                final SAMRecord rec = createSamRecordWithAdapterSequence(readLength, adapterPair, adapterPosition);
                final AdapterPair matchedAdapter = adapterTrimIlluminaSingleRead(rec, minAdapterLength, 0.1, adapterPair);
                final Object xt = rec.getAttribute(XT);

                rec.setAttribute(XT, null);
                final AdapterPair truncatedAdapter = marker.adapterTrimIlluminaSingleRead(rec, minAdapterLength, 0.1);
                final Object xtFromMarker = rec.getAttribute(XT);

                if (adapterPosition <= readLength - minAdapterLength) {
                    assertEquals(matchedAdapter, adapterPair, "Failed at position " + adapterPosition);
                    assertEquals((Integer) xt, valueOf(adapterPosition + 1));
                    assertNotNull(truncatedAdapter);
                    assertEquals((Integer) xtFromMarker, valueOf(adapterPosition + 1));
                } else {
                    assertNull(matchedAdapter);
                    assertNull(xt);
                    assertNull(truncatedAdapter);
                    assertNull(xtFromMarker);
                }
            }
        }
    }

    private SAMRecord createSamRecordWithAdapterSequence(final int readLength, final IlluminaAdapterPair adapterPair, final int adapterPosition) {
        final String adapterString = adapterPair.get3PrimeAdapterInReadOrder();
        final int replacementLength = min(adapterString.length(), readLength - adapterPosition);
        final String adapterSubstring = adapterString.substring(0, replacementLength);
        final String readBases = replaceSubstring(makeBogusReadString(readLength), adapterSubstring, adapterPosition, adapterSubstring.length());
        final SAMRecord rec = new SAMRecord(null);
        rec.setReadString(readBases);
        return rec;
    }

    @DataProvider(name = "testAdapterInAllReadPositionsDataProvider")
    public Object[][] testAdapterInAllReadPositionsDataProvider() {
        return new Object[][]{{100}, {36}};
    }

    @DataProvider(name = "clipTestData")
    public Object[][] getClipTestData() {
        final String FORWARD = PAIRED_END.get3PrimeAdapter();
        final String SE_FORWARD = SINGLE_END.get3PrimeAdapter();
        final String REVERSE = PAIRED_END.get5PrimeAdapterInReadOrder();
        return new Object[][]{
                {"Simple test 1", "AAAAACCCCCAGATCGGAAGAGCA", "AGATCGGAAGAGCG", 14, 0.15, 10},
                {"Simple test 2", "AAAAACCCCCGGGGGAGATCGGAAGAGCA", "AGATCGGAAGAGCG", 14, 0.15, 15},
                {"No adapter", "AAAAACCCCCGGGGGTTTTT", "AGATCGGAAGAGCG", 6, 0.15, -1},
                {"Partial adapter", "AAAAACCCCCTGATCGGAA", "AGATCGGAAGAGCG", 9, 0.15, 10},
// no longer support clips in middle of read
//          {"Adapter+Primer", "AAAAACCCCCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", "AGATCGGAAGAGCG", 6, 0.15, 10},
                {"No sequence", null, "AGATCGGAAGAGCG", 6, 0.15, -1},
                {"Read too short", "AGATCGGAAG", "AGATCGGAAGAGCG", 11, 0.15, -1},
                {"All No Calls", "AAACCCNNNNNNNNNNNNNNNNNNNN", "AGATCGGAAGAGCG", 6, 0.15, -1},
                {"From Test Data1", "CGGCATTCCTGCTGAACCGAGATCGGAAGAGCGTCGTGTAGGGAAAGGGGGTGGATCTCGGTGGGCGGCGTGTTGT", REVERSE, 57, 0.15, 19},
                {"From Test Data2a", "CGGAAGAGCGGTTCAGCAGGAATGCCGAGATCCGAA", REVERSE, 9, 0.14, 27},  // XT:i:28
                {"From Test Data2b", "CGGAAGAGCGGTTCAGCAGGAATGCCGAGATCGGAA", REVERSE, 10, 0.14, -1},  // only matches 9
                {"From PE Test Data1", "CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGGT", FORWARD, 17, 0.14, 19},
                {"From PE Test Data2", "CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGGT", REVERSE, 5, 0.14, -1},
                new Object[]{"From Test 8-clipped", "TGGGGTGGTTATTGTTGATTTTTGTTTGTGTGTTAGGTTGTTTGTGTTAGTTTTTTATTTTATTTTCGAGATCGAA", FORWARD, 8, 0.14, 68},
                {"50% can be bad", "AAAAACCCCCAGGTCGGAAGAGCG", "AGATCGGAAGAGCG", 5, 0.5, 18},        // 18?

                {"From 30E54AAXX.5.a", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", 34, 0.14, 0},  // 0!? From KT test case
                {"From 30E54AAXX.5.b", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", SE_FORWARD, 34, 0.14, 0},  // 1?? From KT test case
        };
    }

    @DataProvider(name = "clipPairedTestData")
    public Object[][] getClipPairedTestData() {
        return new Object[][]{
                {"Basic positive paired test matching", "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "GCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCGG", INDEXED},
                // This matches on first end and matches at higher stringency on that one side
                {"Basic positive first end one-sided test matching", "CTACTGGCGCTGAAAAGATCGGAAGAGCGGTTCAGC", "AAAAAATTTTTTCCCCCCGGGGGGAAAAAATTTTTT", PAIRED_END},
                // This matches on second end and matches at higher stringency on that one side
                {"Basic positive second end one-sided test matching", "AAAAAATTTTTTCCCCCCGGGGGGAAAAAATTTTTT", "AAAAAATTTTTTCCCAGATCGGAAGAGCGTCGTGTA", PAIRED_END},
                // This matches on one side and does not match at higher stringency on that one side
                {"Basic negative one-sided test matching", "GGCGCTGAAACTACTGGCGCTGAAAAGATCGGAAGA", "AAAAAATTTTTTCCCCCCGGGGGGAAAAAATTTTTT", null},
                // These match but at different positions.  No clip should be done.
                {"Mis-match paired test matching", "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "AGCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCG", null},
        };
    }

    private static class SingleEndAdapter implements AdapterPair {
        private final String name;
        private final String threePrimeAdapter;

        private SingleEndAdapter(final String name, final String threePrimeAdapter) {
            this.name = name;
            this.threePrimeAdapter = threePrimeAdapter;
        }

        @Override
        public String get3PrimeAdapter() {
            return threePrimeAdapter;
        }

        @Override
        public String get3PrimeAdapterInReadOrder() {
            return threePrimeAdapter;
        }

        @Override
        public byte[] get3PrimeAdapterBytes() {
            return stringToBytes(threePrimeAdapter);
        }

        @Override
        public byte[] get3PrimeAdapterBytesInReadOrder() {
            return get3PrimeAdapterBytes();
        }

        @Override
        public String get5PrimeAdapter() {
            throw new UnsupportedOperationException();
        }

        @Override
        public String get5PrimeAdapterInReadOrder() {
            throw new UnsupportedOperationException();
        }

        @Override
        public byte[] get5PrimeAdapterBytes() {
            throw new UnsupportedOperationException();
        }

        @Override
        public byte[] get5PrimeAdapterBytesInReadOrder() {
            throw new UnsupportedOperationException();
        }

        @Override
        public String getName() {
            return name;
        }

    }

    /**
     * @return A String of given length with 6 As, 6 Cs, 6 Gs, etc. until length is reached.
     */
    private static String makeBogusReadString(final int len) {
        final StringBuilder builder = new StringBuilder(len);
        final Map<Character, Character> nextChar = new HashMap<>();
        nextChar.put('A', 'C');
        nextChar.put('C', 'G');
        nextChar.put('G', 'T');
        nextChar.put('T', 'A');
        for (char curChar = 'A'; true; curChar = nextChar.get(curChar)) {
            for (int i = 0; i < 6; ++i) {
                if (builder.length() == len) return builder.toString();
                builder.append(curChar);
            }
        }
    }

    private static String patchAdapterSubsequenceIntoRead(final String read, final String adapter, final int length) {
        return replaceSubstring(read, adapter.substring(0, length), read.length() - length, read.length());
    }

    private static String replaceSubstring(final String s, final String replacement, final int position, final int charsToReplace) {
        final StringBuilder sb = new StringBuilder(s);
        sb.replace(position, position + charsToReplace, replacement);
        return sb.toString();
    }


    @Test
    public void testAdapterTruncation() {
        final AdapterMarker marker = new AdapterMarker(30,
                INDEXED,
                DUAL_INDEXED,
                NEXTERA_V2,
                FLUIDIGM);
        // INDEXED and DUAL_INDEXED should collapse at length 30
        assertTrue(marker.getAdapters().length < 4, "Did not collapse: " + marker.getAdapters().length);
    }

    /**
     * Confirm that after the requisite number of sightings of adapter, the list is trimmed
     */
    @Test(dataProvider = "testAdapterListTruncationDataProvider")
    public void testAdapterListTruncation(final IlluminaAdapterPair adapterPair) {
        final int thresholdForTruncatingAdapterList = 20;
        final int readLength = 100;

        // Throw all the adapter pairs into the marker.  There is some danger in doing this because
        // IlluminaAdapterPair.ALTERNATIVE_SINGLE_END has 3' that is a substring of some of the others.  In the enum it appears
        // last so it is tested last.
        final AdapterMarker marker =
                new AdapterMarker(values()).
                        setThresholdForSelectingAdaptersToKeep(thresholdForTruncatingAdapterList);
        int adapterPosition = 1;

        assertTrue(adapterPosition + thresholdForTruncatingAdapterList < readLength - MIN_MATCH_BASES,
                "Test is configured improperly -- eventually will not be enough adapter bases in read.");

        int originalNumberOfAdapters = marker.getAdapters().length;
        for (int i = 0; i < thresholdForTruncatingAdapterList; ++i) {

            // First, a read with no adapter in it.
            SAMRecord rec = new SAMRecord(null);
            rec.setReadString(makeBogusReadString(readLength));
            AdapterPair matchedPair = marker.adapterTrimIlluminaSingleRead(rec);
            assertNull(matchedPair);
            assertNull(rec.getAttribute(XT));

            // Adapter list should not be truncated yet
            assertEquals(marker.getAdapters().length, originalNumberOfAdapters);

            // Then, a record with some adapter sequence in it.
            rec = createSamRecordWithAdapterSequence(readLength, adapterPair, adapterPosition);
            matchedPair = marker.adapterTrimIlluminaSingleRead(rec);
            assertNotNull(matchedPair);
            assertEquals(rec.getIntegerAttribute(XT).intValue(), adapterPosition + 1,
                    rec.getReadString() + " matched " + matchedPair);

            // Put adapter in different place next time.
            adapterPosition++;
        }

        assertEquals(marker.getAdapters().length, 1, "Did not truncate adapter list to 1 element");
        assertTrue(marker.getAdapters()[0].getName().contains(adapterPair.getName()),
                format("Expected '%s' to contain '%s'", marker.getAdapters()[0].getName(), adapterPair.getName()));
    }

    @DataProvider(name = "testAdapterListTruncationDataProvider")
    public Object[][] testAdapterListTruncationDataProvider() {
        Object[][] ret = new Object[values().length][];
        for (int i = 0; i < ret.length; ++i) {
            ret[i] = new Object[]{values()[i]};
        }
        return ret;
    }
}
