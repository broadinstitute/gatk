package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.function.Function;
import java.util.stream.Stream;

public class SplitReadsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DATA_PREFIX = "split_reads";
    private static final String REFERENCE_SEQUENCE = TEST_DATA_PREFIX + ".fasta";

    private File getReferenceSequence() {
        return new File(getTestDataDir(), REFERENCE_SEQUENCE);
    }

    private boolean isReferenceRequired(final SamReader.Type type) {
        return type == SamReader.Type.CRAM_TYPE;
    }

    @DataProvider(name = "splitReadsData", parallel = true)
    public Object[][] getSplitReadsData() {
        final Map<String, Integer> byNone = new TreeMap<>();
        byNone.put("", 19);

        final Map<String, Integer> bySample = new TreeMap<>();
        bySample.put(".Momma", 17);
        bySample.put(".Poppa", 2);

        final Map<String, Integer> byRG = new TreeMap<>();
        byRG.put(".0", 17);
        byRG.put(".1", 2);

        final Map<String, Integer> bySampleAndRG = new TreeMap<>();
        bySampleAndRG.put(".Momma.0", 17);
        bySampleAndRG.put(".Poppa.1", 2);

        final Function<SamReader.Type, Stream<Object[]>> argTests = t -> Stream.of(
                new Object[]{t, Collections.<String>emptyList(), byNone},
                new Object[]{t, Collections.singletonList(SplitReads.SAMPLE_SHORT_NAME), bySample},
                new Object[]{t, Collections.singletonList(SplitReads.READ_GROUP_SHORT_NAME), byRG},
                new Object[]{t, Arrays.asList(SplitReads.SAMPLE_SHORT_NAME, SplitReads.READ_GROUP_SHORT_NAME), bySampleAndRG}
        );

        return getSamReaderTypes()
                .filter(t -> t != SamReader.Type.CRAM_TYPE) // https://github.com/samtools/htsjdk/issues/148 && https://github.com/samtools/htsjdk/issues/153
                .map(argTests)
                .flatMap(Function.identity())
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "splitReadsData")
    public void testSplitReadsByReadGroup(final SamReader.Type type,
                                          final List<String> splitArgs,
                                          final Map<String, Integer> splitCounts) throws Exception {
        final String fileExtension = "." + type.fileExtension();
        final List<String> args = new ArrayList<>();

        Path outputDir = Files.createTempDirectory(
                splitArgs.stream().reduce(TEST_DATA_PREFIX, (acc, arg) -> acc + "." + arg) + fileExtension + "."
        );
        outputDir.toFile().deleteOnExit();

        args.add("-"+ StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(getTestDataDir() + "/" + TEST_DATA_PREFIX + fileExtension);

        args.add("-"+ StandardArgumentDefinitions.OUTPUT_SHORT_NAME );
        args.add(outputDir.toString());

        if (isReferenceRequired(type)) {
            args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME );
            args.add(getTestDataDir()+ "/" + REFERENCE_SEQUENCE);
        }

        splitArgs.forEach(arg -> {
            args.add("-" + arg );
        });

        Assert.assertNull(runCommandLine(args));

        for (final Map.Entry<String, Integer> splitCount: splitCounts.entrySet()) {
            final String outputFileName = TEST_DATA_PREFIX + splitCount.getKey() + fileExtension;
            Assert.assertEquals(
                    getReadCounts(outputDir, outputFileName),
                    (int)splitCount.getValue(),
                    "unexpected read count for " + outputFileName);
        }
    }

    private int getReadCounts(final Path tempDirectory, final String fileName) {
        final File path = tempDirectory.resolve(fileName).toFile();
        IOUtil.assertFileIsReadable(path);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(getReferenceSequence()).open(path);
        int count = 0;
        for (@SuppressWarnings("unused") final SAMRecord rec : in) {
            count++;
        }
        CloserUtil.close(in);
        return count;
    }

    private static Stream<SamReader.Type> getSamReaderTypes() {
        return Stream
                .of(SamReader.Type.class.getFields())
                .filter(f -> Modifier.isStatic(f.getModifiers()))
                .filter(f -> f.getType().isAssignableFrom(SamReader.Type.class))
                .map(f -> orNull(() -> f.get(null)))
                .filter(v -> v instanceof SamReader.Type)
                .map(v -> (SamReader.Type) v);
    }

    private static <V> V orNull(final Callable<V> callable) {
        try {
            return callable.call();
        } catch (Exception e) {
            return null;
        }
    }
}
