package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Compares the VRS IDs GVS computes with those of vrs-python, GA4GH's reference implementation (VS-2059).
 *
 * <p>The expected values are in {@code VrsPythonConcordance/expected.tsv}, one row per ALT, written by
 * {@code generate_expected.py} beside it, which pins the vrs-python version. Each row is checked through
 * {@link VetCreator#computeVrsAlleleRecords}, the path ingest takes, against the checked-in chr20/21 subset
 * of the GVS reference. The normalized start, end and state are asserted before the location and allele IDs,
 * so a failure names the earliest step that differs.</p>
 *
 * <p>A row with a {@code gvs_differs} reason is one where GVS deliberately does not match vrs-python, because
 * the GVS reference N-masks sequence that RefSeq has. Its allele ID is asserted to differ, so a change that
 * makes it match shows up here and the reason can be retired. GVS's own normalized interval and state are
 * asserted too, so the row also guards what GVS returns there.</p>
 */
public class VrsPythonConcordanceUnitTest {
    private static final Path EXPECTED = Paths.get(GATKBaseTest.publicTestDir,
            "org/broadinstitute/hellbender/tools/gvs/ingest/VrsPythonConcordance/expected.tsv");
    private static final Path REFERENCE = Paths.get(GATKBaseTest.largeFileTestDir, "Homo_sapiens_assembly38.20.21.fasta");
    private static final String SPAN_DEL = "*";

    private static final MessageType PARQUET_SCHEMA = MessageTypeParser.parseMessageType("""
            message VariantRow {
            	required int64 sample_id;
            	required int64 location;
            	required binary ref (UTF8);
            	required binary alt (UTF8);
            	required binary call_GT (UTF8);
            	required int64 call_GQ;
            }
            """);

    private ReferenceDataSource reference;
    private VetCreator vetCreator;

    /** One fixture row: a single ALT of a VCF record and vrs-python's result for it. */
    private record Expected(Map<String, String> fields) {
        String get(final String column) { return fields.get(column); }
        String label() { return get("case") + " ALT " + get("alt_index") + " (" + get("alt") + ")"; }
        @Override public String toString() { return get("alt"); }  // keeps TestNG's per-case names short
    }

    @BeforeClass
    public void setUp() throws IOException {
        reference = ReferenceDataSource.of(REFERENCE);
        vetCreator = new VetCreator("vrs_concordance", 100L, "001",
                Files.createTempDirectory("vrs_concordance").toFile(), CommonCode.OutputType.PARQUET,
                "test", "test", true, false, PARQUET_SCHEMA,
                true,   // enableVrsIds
                null);  // vrsAlleleCreator
    }

    @AfterClass(alwaysRun = true)
    public void tearDown() {
        if (vetCreator != null) {
            vetCreator.closeTool();
        }
        if (reference != null) {
            reference.close();
        }
    }

    @DataProvider(name = "cases")
    public Object[][] cases() throws IOException {
        final List<String> lines = Files.readAllLines(EXPECTED).stream().filter(l -> !l.startsWith("#")).toList();
        final List<String> header = Arrays.asList(lines.get(0).split("\t", -1));
        final Map<String, List<Expected>> byCase = new LinkedHashMap<>();
        for (final String line : lines.subList(1, lines.size())) {
            final String[] values = line.split("\t", -1);
            Assert.assertEquals(values.length, header.size(), "malformed fixture row: " + line);
            final Map<String, String> fields = new LinkedHashMap<>();
            for (int i = 0; i < values.length; i++) {
                fields.put(header.get(i), values[i]);
            }
            byCase.computeIfAbsent(fields.get("case"), k -> new ArrayList<>()).add(new Expected(fields));
        }
        return byCase.entrySet().stream()
                .map(e -> new Object[]{e.getKey(), e.getValue()})
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "cases")
    public void testMatchesVrsPython(final String caseName, final List<Expected> alts) {
        final Expected first = alts.get(0);
        final String contig = first.get("contig");
        final int pos = Integer.parseInt(first.get("pos"));
        final String ref = first.get("ref");

        // The fixture's REF must be what the GVS reference has, or the case does not test what it claims.
        final ReferenceContext referenceContext = new ReferenceContext(reference, new SimpleInterval(contig, pos, pos));
        Assert.assertEquals(new String(referenceContext.getBases(new SimpleInterval(contig, pos, pos + ref.length() - 1))),
                ref, caseName + ": REF is not the GVS reference sequence");

        final List<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(ref, true));
        for (final Expected alt : alts) {
            alleles.add(SPAN_DEL.equals(alt.get("alt")) ? Allele.SPAN_DEL : Allele.create(alt.get("alt")));
        }
        alleles.add(Allele.NON_REF_ALLELE);  // as in a gVCF; it gets no entry
        final VariantContext variant = new VariantContextBuilder("vrs_concordance", contig, pos,
                pos + ref.length() - 1, alleles).make();

        final List<VrsAlleleRecord> records = vetCreator.computeVrsAlleleRecords(variant, referenceContext);
        final List<String> vrsAlleleIds = VetCreator.vrsAlleleIds(records);
        Assert.assertEquals(records.size(), alts.size(), caseName + ": one record per ALT");

        for (final Expected expected : alts) {
            final int index = Integer.parseInt(expected.get("alt_index"));
            final String label = expected.label();
            final VrsAlleleRecord actual = records.get(index);

            if (SPAN_DEL.equals(expected.get("alt"))) {
                Assert.assertNull(actual, label + ": * has no VRS representation");
                Assert.assertEquals(vrsAlleleIds.get(index), VetCreator.NO_VRS_ALLELE_ID, label);
                continue;
            }
            Assert.assertNotNull(actual, label);

            if (!expected.get("gvs_differs").isEmpty()) {
                Assert.assertNotEquals(actual.vrsAlleleId, expected.get("allele_id"),
                        label + " now matches vrs-python; retire its gvs_differs reason: " + expected.get("gvs_differs"));
                // Not matching vrs-python is not enough: pin what GVS itself returns. A substitution over N
                // has nothing to trim or roll, so GVS keeps the input interval and writes the ALT literally.
                Assert.assertEquals(expected.get("alt").length(), ref.length(),
                        label + ": only substitutions have a GVS expectation here; add one for this kind of row");
                Assert.assertEquals(actual.start, pos - 1L, label + ": GVS start");
                Assert.assertEquals(actual.end, pos - 1L + ref.length(), label + ": GVS end");
                Assert.assertEquals(actual.stateType, VrsAlleleRecord.STATE_LITERAL, label + ": GVS state type");
                Assert.assertEquals(actual.stateSequence, expected.get("alt"), label + ": GVS state sequence");
                Assert.assertEquals(vrsAlleleIds.get(index), actual.vrsAlleleId, label + ": vrs_allele_ids entry");
                continue;
            }

            // Normalization first, then serialization, so the first failure names the step that differs.
            Assert.assertEquals(actual.start, Long.parseLong(expected.get("start")), label + ": normalized start");
            Assert.assertEquals(actual.end, Long.parseLong(expected.get("end")), label + ": normalized end");
            Assert.assertEquals(actual.stateType, expected.get("state_type"), label + ": state type");
            if (VrsAlleleRecord.STATE_LITERAL.equals(actual.stateType)) {
                Assert.assertEquals(actual.stateSequence, expected.get("state_sequence"), label + ": state sequence");
            } else {
                Assert.assertEquals(actual.stateLength, Integer.valueOf(expected.get("state_length")),
                        label + ": state length");
                Assert.assertEquals(actual.stateRepeatSubunitLength,
                        Integer.valueOf(expected.get("state_repeat_subunit_length")),
                        label + ": state repeatSubunitLength");
            }
            Assert.assertEquals(actual.vrsLocationId, expected.get("location_id"), label + ": location ID");
            Assert.assertEquals(actual.vrsAlleleId, expected.get("allele_id"), label + ": allele ID");
            Assert.assertEquals(vrsAlleleIds.get(index), expected.get("allele_id"), label + ": vrs_allele_ids entry");
        }
    }
}
