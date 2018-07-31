package org.broadinstitute.hellbender.utils.codecs.table;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static java.util.Arrays.asList;
import static java.util.Collections.emptyList;

public final class TableCodecUnitTest extends GATKBaseTest {

    @DataProvider(name = "badNames")
    public Object[][] badNames() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{"a.tsv"});
        params.add(new String[]{"a.table.gz"});
        params.add(new String[]{"a.bed"});
        params.add(new String[]{"a.bcf"});
        params.add(new String[]{"a.hapmap"});
        params.add(new String[]{"a.refseq"});
        params.add(new String[]{"a.beagle"});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "badNames")
    public void testBadNames(String badName){
        TableCodec tc = new TableCodec();
        Assert.assertFalse(tc.canDecode(badName), badName);
    }

    @DataProvider(name = "goodNames")
    public Object[][] goodNames() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{"a.table"});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "goodNames")
    public void testGoodNames(String goodName){
        TableCodec tc = new TableCodec();
        Assert.assertTrue(tc.canDecode(goodName), goodName);
    }

    @Test
    public void testChrs(){
        TableCodec tc = new TableCodec();
        Assert.assertEquals(tc.decode("1:1  1   2   3").getContig(), "1");
        Assert.assertEquals(tc.decode("chr1:1  1   2   3").getContig(), "chr1");
        Assert.assertEquals(tc.decode("1:1+  1   2   3").getContig(), "1");
        Assert.assertEquals(tc.decode("1  1   2   3").getContig(), "1");
        Assert.assertEquals(tc.decode("fred  1   2   3").getContig(), "fred");
        Assert.assertEquals(tc.decode("2:1,000  1   2   3").getContig(), "2");
    }

    @DataProvider(name = "dataLines")
    public Object[][] dataLines() {
        List<Object[]> params = new ArrayList<>();
        params.add(new Object[]{"1:1  1   2   3", "1", 1, 1, emptyList(), asList("1:1", "1", "2", "3")});
        params.add(new Object[]{"1:1-2  1   2   3", "1", 1, 2, emptyList(), asList("1:1-2", "1", "2", "3")});
        params.add(new Object[]{"1  1   2   3", "1", 1, Integer.MAX_VALUE, emptyList(), asList("1", "1", "2", "3")});
        params.add(new Object[]{"1:1,000-2,000  1   2   3", "1", 1000, 2000, emptyList(), asList("1:1,000-2,000", "1", "2", "3")});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "dataLines")
    public void testDecode1(String str, String contig, int start, int end, List<String> header, List<String> vals){
        TableCodec tc = new TableCodec();
        final TableFeature decode = tc.decode(str);
        Assert.assertEquals(decode.getContig(), contig, "contig");
        Assert.assertEquals(decode.getStart(), start, "start");
        Assert.assertEquals(decode.getEnd(), end, "end");
        Assert.assertEquals(decode.getHeader(), header, decode.getHeader().toString());

        Assert.assertEquals(decode.columnCount(), vals.size());
        for(int i = 0; i < decode.columnCount(); i++){
            Assert.assertEquals(decode.getValue(i), vals.get(i), "i:" + i);
        }
    }

    @Test
    public void testDecodeHeader(){
        TableCodec tc = new TableCodec();
        LineReader reader= makeReader(asList("HEADER a b c"));
        LineIterator li= new LineIteratorImpl(reader);
        List<String> hd = tc.readActualHeader(li);
        Assert.assertEquals(hd, asList("HEADER", "a", "b", "c"));
    }

    @Test
    public void testDecodeHeader2(){
        TableCodec tc = new TableCodec();
        final String str2= "1:1  1   2   3";
        LineReader reader= makeReader(asList("HEADER a b c", str2));
        LineIterator li= new LineIteratorImpl(reader);
        List<String> hd = tc.readActualHeader(li);
        Assert.assertEquals(hd, asList("HEADER", "a", "b", "c"));

        final TableFeature decode = tc.decode(str2);
        Assert.assertEquals(decode.get("a"), "1");
        Assert.assertEquals(decode.get("b"), "2");
        Assert.assertEquals(decode.get("c"), "3");
        Assert.assertEquals(decode.getLocation().getContig(), "1");
        Assert.assertEquals(decode.getContig(), "1");
        Assert.assertEquals(decode.getLocation().getStart(), 1);
        Assert.assertEquals(decode.getLocation().getEnd(), 1);
    }

    @Test(expectedExceptions = UserException.MalformedFile.class)
    public void testDecodeFailsNoHeader(){
        TableCodec tc = new TableCodec();
        LineReader reader= makeReader(asList("1:1  1   2   3"));
        LineIterator li= new LineIteratorImpl(reader);
        tc.readActualHeader(li);
    }

    @Test
    public void testDecodeOnlyComments(){
        TableCodec tc = new TableCodec();
        LineReader reader= makeReader(asList("#HEADER a b c", "#HEADER d e f"));
        LineIterator li= new LineIteratorImpl(reader);
        final List<String> strings = tc.readActualHeader(li);
        Assert.assertEquals(strings, emptyList());
    }

    @Test
    public void testTwoHeaders(){
        TableCodec tc = new TableCodec();
        LineReader reader= makeReader(asList("HEADER a b c", "HEADER d e f"));
        LineIterator li= new LineIteratorImpl(reader);
        final List<String> strings = tc.readActualHeader(li);
        Assert.assertEquals(strings, asList("HEADER", "a", "b", "c"));
    }

    @Test(expectedExceptions =  UserException.MalformedFile.class)
    public void testTwoHeadersFailsOnRepeat(){
        TableCodec tc = new TableCodec();
        Assert.assertEquals(tc.readActualHeader(new LineIteratorImpl(makeReader(asList("HEADER a b c")))), asList("HEADER", "a", "b", "c"));

        Assert.assertEquals(tc.readActualHeader(new LineIteratorImpl(makeReader(asList("HEADER a b c")))), asList("HEADER", "a", "b", "c"));
    }

    @Test
    public void testDecodeComment(){
        TableCodec tc = new TableCodec();
        LineReader reader= makeReader(asList("#HEADER a b c", "HEADER d e f"));
        LineIterator li= new LineIteratorImpl(reader);
        List<String> hd = tc.readActualHeader(li);
        Assert.assertEquals(hd, asList("HEADER", "d", "e", "f"));
    }

    private LineReader makeReader(List<String> strings) {
        return new LineReader() {
            private Iterator<String> iterator = strings.iterator();
            @Override
            public String readLine() throws IOException {
                return iterator.hasNext() ? iterator.next() : null;
            }

            @Override
            public void close() {
            }
        };
    }

    @DataProvider(name = "stringNull")
    public Object[][] stringNull() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{TableCodec.HEADER_DELIMITER + " foo"});
        params.add(new String[]{TableCodec.HEADER_DELIMITER + " bar"});
        params.add(new String[]{TableCodec.IGV_HEADER_DELIMITER + " baz"});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "stringNull")
    public void testDecodeNull(String stringNull){
        TableCodec tc = new TableCodec();
        final TableFeature decode = tc.decode(stringNull);
        Assert.assertNull(decode, stringNull);
    }

    @DataProvider(name = "stringNotNull")
    public Object[][] stringNotNull() {
        List<Object[]> params = new ArrayList<>();
        params.add(new String[]{"foo " + TableCodec.HEADER_DELIMITER});
        params.add(new String[]{"foo " + TableCodec.HEADER_DELIMITER});
        params.add(new String[]{"foo " + TableCodec.IGV_HEADER_DELIMITER});
        return params.toArray(new Object[][]{});
    }

    @Test(dataProvider = "stringNotNull")
    public void testDecodeNotNull(String stringNotNull){
        TableCodec tc = new TableCodec();
        final TableFeature decode = tc.decode(stringNotNull);
        Assert.assertNotNull(decode, stringNotNull);
    }
}
