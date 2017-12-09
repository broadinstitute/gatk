package org.broadinstitute.hellbender.utils.help;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link GATKHelpDoclet} and {@link GATKHelpDocWorkUnitHandler}.
 */
public class GATKHelpDocletTest extends BaseTest {

    @Test(dataProvider = "testPicardTransformationData")
    public void testPicardTransformationBlock(final String input, final String expected) {
        final StringBuilder original = new StringBuilder(input);
        final GATKHelpDocWorkUnitHandler handler = new GATKHelpDocWorkUnitHandler(new GATKHelpDoclet());
        final StringBuilder result = new StringBuilder(original.length() << 1);
        handler.translatePicardCodeBlock(original, 0, original.length(), result);
        Assert.assertEquals(result.toString(), expected);
    }
    
    @Test
    public void testPicardTransformationAllBlocks() {
    		final String[][] allData = Arrays.stream(testPicardTransformationData()).map(o -> new String[] { "" + o[0], "" + o[1] })
    				.toArray(String[][]::new);
    		final GATKHelpDocWorkUnitHandler handler = new GATKHelpDocWorkUnitHandler(new GATKHelpDoclet());
    		final String input = Arrays.stream(allData).map(s -> s[0]).collect(Collectors.joining("In the middle", "Start", "End"));
    		final String expected = Arrays.stream(allData).map(s -> s[1]).collect(Collectors.joining("In the middle", "Start", "End"));
    		Assert.assertEquals(handler.translatePicardCodeBlocks(input), expected);
    }
    
    @Test(dataProvider = "testPicardTransformationData")
    public void testPicardTransformationBlocks(final String input, final String expected) {
        final StringBuilder original = new StringBuilder(input);
        final GATKHelpDocWorkUnitHandler handler = new GATKHelpDocWorkUnitHandler(new GATKHelpDoclet());
        final String result = handler.translatePicardCodeBlocks(original);
        Assert.assertEquals(result, expected);
    }

    @DataProvider(name = "testPicardTransformationData")
    public Object[][] testPicardTransformationData() {

        final String[] interleavedInputExceptedPairs = {
                // simple one liner with a flag.
                lines("<pre>java -jar picard.jar MergeVcfs I=file1 O=file2 CREATE_INDEX=true INPUT=file3</pre>"),
                lines("<pre>gatk MergeVcfs -I file1 -O file2 --CREATE_INDEX --INPUT file3</pre>"),
                // same in multiple lines:
                lines("<pre>   ",
                      "        java -jar picard.jar MergeVcfs \\  ",
                      "                             I=file1   \\",
                      "                             O=file2   \\",
                      "                             INPUT=fil3 \\",
                      "                             CREATE_INDEX=true  ",
                      "</pre>"),
                lines("<pre>   ",
                        "        gatk MergeVcfs \\  ",
                        "                             -I file1   \\",
                        "                             -O file2   \\",
                        "                             --INPUT fil3 \\",
                        "                             --CREATE_INDEX  ",
                        "</pre>"),
                // <pre> share line with other elements of the command line example
                lines("<pre>   java -jar picard.jar MergeVcfs \\  ",
                        "                             I=file1   \\",
                        "                             O=file2   \\",
                        "                             INPUT=fil3 \\",
                        "                             CREATE_INDEX=true </pre>"),
                lines("<pre>   gatk MergeVcfs \\  ",
                        "                             -I file1   \\",
                        "                             -O file2   \\",
                        "                             --INPUT fil3 \\",
                        "                             --CREATE_INDEX </pre>"),
                // mulitple commands:
                lines("<pre>   java -jar picard.jar MergeVcfs \\  ",
                        "                             I=file1  ",
                        "      java -jar picard.jar MergeVcfs O=file2   \\",
                        "                             INPUT=fil3 \\",
                        "                             CREATE_INDEX=true </pre>"),
                lines("<pre>   gatk MergeVcfs \\  ",
                        "                             -I file1  ",
                        "      gatk MergeVcfs -O file2   \\",
                        "                             --INPUT fil3 \\",
                        "                             --CREATE_INDEX </pre>"),
                lines("<pre>",
                        "java -jar picard.jar FastqToSam \\",
                        "    F1=input_reads.fastq \\",
                        "    O=unaligned_re ads.bam \\",
                        "    SM=sample001 \\",
                        "    RG=rg0013",
                        "</pre>"),
                lines("<pre>",
                        "gatk FastqToSam \\",
                        "    -F1 input_reads.fastq \\",
                        "    -O unaligned_re ads.bam \\",
                        "    -SM sample001 \\",
                        "    -RG rg0013",
                        "</pre>"),
                lines("<pre>",
                        "java -Xmx10g -Djava.io.tmp=mytemp   -jar picard.jar FastqToSam \\",
                        "    F1=input_reads.fastq \\",
                        "    O=unaligned_re ads.bam \\",
                        "    SM=sample001 \\",
                        "    RG=rg0013",
                        "</pre>"),
                lines("<pre>",
                        "gatk --javaOptions '-Xmx10g -Djava.io.tmp=mytemp' FastqToSam \\",
                        "    -F1 input_reads.fastq \\",
                        "    -O unaligned_re ads.bam \\",
                        "    -SM sample001 \\",
                        "    -RG rg0013",
                        "</pre>"),
                lines("<pre>",
                		"some other commands",
                		"some other commands 2",
                		"java -jar picard.jar FastqToSam \\",
                		"     F1=input_reads.fastq \\",
                		"     O=output.bam",
                		"some other commands 3",
                		"",
                		"",
                		"</pre>"),
                lines("<pre>",
                		"some other commands",
                		"some other commands 2",
                		"gatk FastqToSam \\",
                		"     -F1 input_reads.fastq \\",
                		"     -O output.bam",
                		"some other commands 3",
                		"",
                		"",
                		"</pre>"),
                lines("sddasd<pre>sddddsd",
                		"java -jar \t random-picard-0.0.1_v13.jar MergeVcfs</pre>"),
                lines("sddasd<pre>sddddsd",
                		"gatk MergeVcfs</pre>"),
                lines("<pre > java -jar picard.jar MergeVcfs </pre >"),
                lines("<pre > gatk MergeVcfs </pre >")
                
        };

        final List<Pair<String, String>> result = new ArrayList<>();
        for (int i = 0; i < interleavedInputExceptedPairs.length; i += 2) {
            result.add(new ImmutablePair<>(interleavedInputExceptedPairs[i], interleavedInputExceptedPairs[i + 1]));
        }
        return result.stream().map(pair -> new Object[] { pair.getLeft(), pair.getRight()}).toArray(Object[][]::new);
    }

    private static String lines(String ... lines) {
        return Arrays.stream(lines).collect(Collectors.joining("\n"));
    }
}
