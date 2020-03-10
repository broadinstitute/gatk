package org.broadinstitute.hellbender.utils.codecs.gtf;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Test class for the {@link AbstractGtfCodec}.
 * Modeled after the {@link org.broadinstitute.hellbender.utils.codecs.table.TableCodecUnitTest}, with extras specific to this file format.
 * Created by jonn on 7/27/17.
 */
public class AbstractGtfCodecUnitTest extends GATKBaseTest {

    List<String> prependCommentToArrayElements(final String[] elements, final String comment) {
        return Arrays.stream(elements).map(s -> comment + s).collect(Collectors.toList());
    }

    @DataProvider
    Object[][] provideForCheckHeaderLineStartsWith() {

        final GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();

        return new Object[][] {
                {gencodeGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, gencodeGtfCodec.getDefaultLineComment()), 0, "LINE_1", true},
                {gencodeGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, gencodeGtfCodec.getDefaultLineComment()), 0, "LINE_2", false},
                {ensemblGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, ensemblGtfCodec.getDefaultLineComment()), 0, "LINE_1", true},
                {ensemblGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, ensemblGtfCodec.getDefaultLineComment()), 0, "LINE_2", false},
        };
    }

    @Test(dataProvider = "provideForCheckHeaderLineStartsWith")
    void testCheckHeaderLineStartsWith(final AbstractGtfCodec codec, final List<String> header, final int lineNum,
                                       final String startingText, final boolean expected) {
        Assert.assertEquals(codec.checkHeaderLineStartsWith(header, lineNum, startingText), expected);
    }

    @DataProvider
    Object[][] provideForCheckHeaderLineStartsWith_WithThrow() {
        final GencodeGtfCodec gencodeGtfCodec = new GencodeGtfCodec();
        final EnsemblGtfCodec ensemblGtfCodec = new EnsemblGtfCodec();

        return new Object[][] {
                {gencodeGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, gencodeGtfCodec.getDefaultLineComment()), 0, "LINE_2"},
                {ensemblGtfCodec, prependCommentToArrayElements(new String[] {"LINE_1", "LINE_2", "LINE_3", "LINE_4", "LINE_5"}, ensemblGtfCodec.getDefaultLineComment()), 0, "LINE_2"},
        };
    }

    @Test(dataProvider = "provideForCheckHeaderLineStartsWith_WithThrow", expectedExceptions = UserException.MalformedFile.class)
    void testCheckHeaderLineStartsWith_WithThrow(final AbstractGtfCodec codec, final List<String> header, final int lineNum,
                                       final String startingText) {
        codec.checkHeaderLineStartsWith(header, lineNum, startingText, true);
    }
}
