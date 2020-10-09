package org.broadinstitute.hellbender.tools.dragstr;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorTestConstants;
import org.broadinstitute.hellbender.utils.BinaryTableReader;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dragstr.STRTableFile;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;


/**
 * Integration tests for {@link ComposeSTRTableFile}.
 *
 * <p>
 *     What we did here is to use DRAGEN 3.4.? boxes to generate the corresponding str_table.bin file for
 *     two reference present in src/test/resources, transcribe it into a easier to read in tsv table and compare it
 *     with what we get
 * </p>
 */
public final class ComposeSTRTableFileIntegrationTest extends CommandLineProgramTest {

    private static final File EXPECTED_HUMAN_B37_20_21_TBL_GZ =
            new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/compose-str-table-file-human_g1k_v37.20.21.dragen.tab.gz");
    private static final File EXPECTED_HUMAN_B38_20_21_TBL_GZ =
            new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/compose-str-table-file-human_v38.20.21.dragen.tab.gz");
    private static final File EXPECTED_GENCODE_TBL_GZ =
            new File("src/test/resources/large/org/broadinstitute/hellbender/tools/dragstr/compose-str-table-file-gencode-transcripts.dragen.tab.gz");


    private static final File EXPECTED_HUMAN_B37_20_21_REF_DICT =
            new File(b37_reference_20_21.replace(".fasta", ".dict"));

    @Test(dataProvider="testCases")
    public void testTestCase(final TestCase testCase) throws IOException  {
        final File output = File.createTempFile("cstf-test-output", ".zip");
        output.deleteOnExit();
        try {
            final ArgumentsBuilder args = new ArgumentsBuilder();
            args.add(StandardArgumentDefinitions.REFERENCE_SHORT_NAME, testCase.fasta);
            args.addOutput(output);
            runCommandLine(args);
            assertMatchesExpected(output, testCase);
        } finally {
            output.delete();
        }
    }

    @DataProvider
    public Iterator<Object[]> testCases() {
        final List<TestCase> result = new ArrayList<>(3);
        result.add(new TestCase("b37", EXPECTED_HUMAN_B37_20_21_TBL_GZ, new File(b37_reference_20_21)));
        result.add(new TestCase("b38", EXPECTED_HUMAN_B38_20_21_TBL_GZ, new File(b38_reference_20_21)));
        result.add(new TestCase("gencode", EXPECTED_GENCODE_TBL_GZ, new File(FuncotatorTestConstants.GENCODE_TRANSCRIPT_FASTA_FILE_NAME)));
        return result.stream()
                .map(tc -> new Object[] {tc})
                .iterator();
    }

    private void assertMatchesExpected(final File actualFile, final TestCase testCase) throws IOException {
        final STRTableFile actualTable = STRTableFile.open(new GATKPath(actualFile.toString()));
        final SAMSequenceDictionary actualDictionary = actualTable.dictionary();
        final SAMSequenceDictionary expectedDictionary = ReferenceUtils.loadFastaDictionary(testCase.dict);
        Assert.assertSame(SequenceDictionaryUtils.compareDictionaries(actualDictionary, expectedDictionary, true),
                SequenceDictionaryUtils.SequenceDictionaryCompatibility.IDENTICAL);
        final Map<Integer, List<ExpectedLocus>> allExpected;
        try (final ExpectedTableReader reader = new ExpectedTableReader(testCase.dragenOutput)) {
            allExpected = reader.stream()
                  .collect(Collectors.groupingBy(el -> el.chrIdx));
        }
        for (final SAMSequenceRecord seq : actualDictionary.getSequences()) {
            final SimpleInterval si = new SimpleInterval(seq.getSequenceName(), 1, seq.getSequenceLength());
            try (final BinaryTableReader<DragstrLocus> lociReader = actualTable.locusReader(si)) {
                final List<ExpectedLocus> actual = lociReader.stream()
                        .map(ExpectedLocus::new)
                        .collect(Collectors.toList());
                final List<ExpectedLocus> expectedLociInChr = allExpected.getOrDefault(seq.getSequenceIndex(), Collections.emptyList());
                Assert.assertEquals(actual, expectedLociInChr);
            }
        }
    }

    private static class TestCase {
        final String name;
        final File fasta;
        final File dict;
        final File dragenOutput;

        TestCase (final String name, final File expected, final File fasta) {
            this.name = name;
            this.fasta = fasta;
            this.dragenOutput = expected;
            this.dict = new File(fasta.toString().replace(".fasta", ".dict"));
        }

        public String toString() {
            return name;
        }
    }

    private static class ExpectedLocus {
        final int chrIdx;
        final long position;
        final int period;
        final int length;
        final int repeatLength;
        final long mask;

        ExpectedLocus(final DataLine dl) {
            chrIdx = dl.getInt("chrIdx");
            position = dl.getInt("pos");
            period = dl.getInt("period");
            length = dl.getInt("length");
            repeatLength = dl.getInt("repeatLen");
            mask = dl.getInt("mask");
        }

        ExpectedLocus(final DragstrLocus dragstrLocus) {
            chrIdx = dragstrLocus.getChromosomeIndex();
            position = dragstrLocus.getStart();
            period = dragstrLocus.getPeriod();
            repeatLength = Math.min(20, dragstrLocus.getRepeats());
            length = (int) (dragstrLocus.getEnd() - dragstrLocus.getStart() + 1);
            mask = dragstrLocus.getMask();
        }

        @Override
        public int hashCode() {
            return (int) (((((chrIdx * 31) +  position) * 31 + repeatLength) * 31 + length) * 31 + mask) * 31;
        }

        @Override
        public boolean equals(final Object other) {
            if (other instanceof ExpectedLocus) {
                final ExpectedLocus otherCasted = ExpectedLocus.class.cast(other);
                return otherCasted.chrIdx == this.chrIdx && otherCasted.period == this.period && otherCasted.mask == this.mask &&
                        otherCasted.repeatLength == this.repeatLength && otherCasted.length == this.length;
            } else {
                return false;
            }
        }
    }

    private static class ExpectedTableReader extends TableReader<ExpectedLocus> {

        public ExpectedTableReader(File file) throws IOException {
            super(file.toPath());
        }

        @Override
        protected ExpectedLocus createRecord(DataLine dataLine) {
            return new ExpectedLocus(dataLine);
        }
    }
}
