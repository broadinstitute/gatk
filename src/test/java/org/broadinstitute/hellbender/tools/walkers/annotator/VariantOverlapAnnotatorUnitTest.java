package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class VariantOverlapAnnotatorUnitTest extends GATKBaseTest {

    private VariantContext makeVC(final String source, final String id, final List<String> alleles) {
        final VariantContext vc = GATKVariantContextUtils.makeFromAlleles(source, "20", 10, alleles);
        return new VariantContextBuilder(vc).id(id).make();
    }

    private VariantOverlapAnnotator makeAnnotator(final File file, final String dbSNP, final String... overlaps) {
        final FeatureInput<VariantContext> dbSNPBinding = dbSNP == null ? null : new FeatureInput<>(file.getAbsolutePath(), dbSNP, Collections.emptyMap());
        final Map<FeatureInput<VariantContext>, String> overlapBinding = new LinkedHashMap<>();
        for ( final String overlap : overlaps ) {
            overlapBinding.put(new FeatureInput<>(file.getAbsolutePath(), overlap, Collections.emptyMap()), overlap);
        }
        if (overlapBinding.isEmpty()) {                     //to test both constructors
            return new VariantOverlapAnnotator(dbSNPBinding);
        } else {
            return new VariantOverlapAnnotator(dbSNPBinding, overlapBinding);
        }
    }

    @Test
    public void testCreateWithSpecialNames() {
        final String decoyPath= "fred";
        final List<String> names = Arrays.asList("X", "Y", "Z");
        final Map<FeatureInput<VariantContext>, String> overlapBinding = new LinkedHashMap<>();
        for ( final String overlap : names ) {
            overlapBinding.put(new FeatureInput<>(decoyPath, overlap + "Binding", Collections.emptyMap()), overlap);
        }
        final VariantOverlapAnnotator annotator = new VariantOverlapAnnotator(null, overlapBinding);
        Assert.assertEquals(annotator.getOverlapNames(), names);
    }

    @DataProvider(name = "AnnotateRsIDData")
    public Object[][] makeAnnotateRsIDData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        final VariantContext callNoIDAC = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("A", "C"));
        final VariantContext callNoIDAT = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("A", "T"));
        final VariantContext callIDAC = makeVC("call", "foo", Arrays.asList("A", "C"));
        final VariantContext callExistingIDAC = makeVC("call", "rsID1", Arrays.asList("A", "C"));

        final VariantContext dbSNP_AC = makeVC("DBSNP", "rsID1", Arrays.asList("A", "C"));
        final VariantContext dbSNP_AT = makeVC("DBSNP", "rsID2", Arrays.asList("A", "T"));
        final VariantContext dbSNP_AG = makeVC("DBSNP", "rsID3", Arrays.asList("A", "G"));
        final VariantContext dbSNP_AC_AT = makeVC("DBSNP", "rsID1;rsID2", Arrays.asList("A", "C", "T"));
        final VariantContext dbSNP_AC_AG = makeVC("DBSNP", "rsID1;rsID3", Arrays.asList("A", "C", "G"));

        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AC), dbSNP_AC.getID(), true});
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AT), VCFConstants.EMPTY_ID_FIELD, false});
        tests.add(new Object[]{callIDAC, Arrays.asList(dbSNP_AC), "foo" + ";" + dbSNP_AC.getID(), true});
        tests.add(new Object[]{callIDAC, Arrays.asList(dbSNP_AT), "foo", false});
        tests.add(new Object[]{callExistingIDAC, Arrays.asList(dbSNP_AC), "rsID1", true});
        tests.add(new Object[]{callExistingIDAC, Arrays.asList(dbSNP_AT), "rsID1", false});

        final VariantContext callNoIDACT = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("A", "C", "T"));
        tests.add(new Object[]{callNoIDACT, Arrays.asList(dbSNP_AC), dbSNP_AC.getID(), true});
        tests.add(new Object[]{callNoIDACT, Arrays.asList(dbSNP_AT), dbSNP_AT.getID(), true});
        tests.add(new Object[]{callNoIDACT, Arrays.asList(dbSNP_AG), VCFConstants.EMPTY_ID_FIELD, false});
        tests.add(new Object[]{callNoIDACT, Arrays.asList(dbSNP_AC_AT), dbSNP_AC_AT.getID(), true});
        tests.add(new Object[]{callNoIDACT, Arrays.asList(dbSNP_AC_AG), dbSNP_AC_AG.getID(), true});

        // multiple options
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AC, dbSNP_AT), "rsID1", true});
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AT, dbSNP_AC), "rsID1", true});
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AC_AT), "rsID1;rsID2", true});
        tests.add(new Object[]{callNoIDAT, Arrays.asList(dbSNP_AC_AT), "rsID1;rsID2", true});
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AC_AG), "rsID1;rsID3", true});
        tests.add(new Object[]{callNoIDAT, Arrays.asList(dbSNP_AC_AG), VCFConstants.EMPTY_ID_FIELD, false});

        //multiple rsid which match multiallelic site
        final VariantContext callNoIDT_C_TAC = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("T", "C", "TAC"));

        final VariantContext dbSNP_T_TAC_TATAC = makeVC("DBSNP", "rsID1", Arrays.asList("T", "TAC", "TATAC"));
        final VariantContext dbSNP_T_C = makeVC("DBSNP", "rsID2", Arrays.asList("T", "C"));

        tests.add(new Object[]{callNoIDT_C_TAC, Arrays.asList(dbSNP_T_TAC_TATAC, dbSNP_T_C), "rsID1;rsID2", true});
        tests.add(new Object[]{callNoIDT_C_TAC, Arrays.asList(dbSNP_T_C,dbSNP_T_TAC_TATAC), "rsID2;rsID1", true});

        //mixed multiallelic in dbsnp
        final VariantContext callNOID_T_TTCC = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("T", "TTCC"));

        final VariantContext dbSNP_complex_mixed_site = makeVC("DBSNP", "rsID1", Arrays.asList("TTCCTCCTCCTCCTCCTCC", "T", "TTCCTCCTCCTCCTCCTCCTCC"));

        tests.add(new Object[]{callNOID_T_TTCC, Arrays.asList(dbSNP_complex_mixed_site), "rsID1", true});

        //mixed multialleleic in call
        final VariantContext callNOID_complex_mixed_site = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("TTCCTCCTCCTCCTCCTCC", "T", "TTCCTCCTCCTCCTCCTCCTCC"));

        final VariantContext dbSNP_T_TTCC = makeVC("DBSNP", "rsID1", Arrays.asList("T", "TTCC"));

        tests.add(new Object[]{callNOID_complex_mixed_site, Arrays.asList(dbSNP_T_TTCC), "rsID1", true});

        //dbsnp and call both deletions different length
        final VariantContext call_deletion = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("TTTT", "T"));

        final VariantContext dbSNP_deletion_different_length = makeVC("DBSNP", "rsID1", Arrays.asList("TTTTTT", "T"));
        final VariantContext dbSNP_deletion_and_insertion_multiallelic = makeVC("DBSNP", "rsID2", Arrays.asList("TTTTTT", "T", "TTTTTTTTTTT"));

        tests.add(new Object[]{call_deletion, Arrays.asList(dbSNP_deletion_different_length, dbSNP_deletion_and_insertion_multiallelic), VCFConstants.EMPTY_ID_FIELD, false});

        //dbsnp and call both insert different lengths
        final VariantContext call_insertion = makeVC("call", VCFConstants.EMPTY_ID_FIELD, Arrays.asList("T", "TTTT"));

        final VariantContext dbSNP_insertion_different_length = makeVC("DBSNP", "rsID1", Arrays.asList("T", "TTTTTT"));

        tests.add(new Object[]{call_insertion, Arrays.asList(dbSNP_insertion_different_length, dbSNP_deletion_and_insertion_multiallelic), VCFConstants.EMPTY_ID_FIELD, false});

        final VariantContext dbSNP_AC_FAIL = new VariantContextBuilder(makeVC("DBSNP", "rsID1", Arrays.asList("A", "C"))).filter("FAIL").make();
        tests.add(new Object[]{callNoIDAC, Arrays.asList(dbSNP_AC_FAIL), VCFConstants.EMPTY_ID_FIELD, false});


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testAnnotateRsID(final VariantContext toAnnotate, final List<VariantContext> dbSNPRecords, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File("fred");
        final VariantOverlapAnnotator annotator = makeAnnotator(file, "dbsnp");
        final VariantContext annotated = VariantOverlapAnnotator.annotateRsID(dbSNPRecords, toAnnotate);
        Assert.assertNotNull(annotated);
        Assert.assertEquals(annotated.getID(), expectedID);
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testAnnotateOverlaps(final VariantContext toAnnotate, final List<VariantContext> records, final String expectedID, final boolean expectOverlap) throws Exception {
        final String name = "binding";
        final File file= new File("fred");
        final VariantOverlapAnnotator annotator = makeAnnotator(file, null, name);
        final VariantContext annotated = annotator.annotateOverlap(records, name, toAnnotate);
        Assert.assertNotNull(annotated);
        Assert.assertEquals(annotated.getID(), toAnnotate.getID(), "Shouldn't modify annotation");
        Assert.assertEquals(annotated.hasAttribute(name), expectOverlap, "Attribute:" + name);
        if ( expectOverlap ) {
            Assert.assertEquals(annotated.getAttribute(name), true);
        }
    }


    @Test(dataProvider = "AnnotateRsIDData")
    public void testAnnotateOverlapsEmpty(final VariantContext toAnnotate, final List<VariantContext> records, final String expectedID, final boolean expectOverlap) throws Exception {
        final String name = "binding";
        final File file= new File("fred");
        final VariantOverlapAnnotator annotator = makeAnnotator(file, null);
        final VariantContext annotated = annotator.annotateOverlap(records, name, toAnnotate);
        Assert.assertEquals(annotated, toAnnotate);
    }

    @CommandLineProgramProperties(summary = "", oneLineSummary = "", programGroup=TestProgramGroup.class)
    private static class ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap extends CommandLineProgram {
        @Argument(fullName = "dbsnp", shortName = "f")
        FeatureInput<Feature> featureArgument;

        @Argument(fullName = "binding", shortName = "b")
        FeatureInput<Feature> binding;

        public ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(File f) {
            featureArgument = new FeatureInput<>(f.getAbsolutePath(), "dbsnp", Collections.emptyMap());
            binding = new FeatureInput<>(f.getAbsolutePath(), "binding", Collections.emptyMap());
        }

        @Override
        protected Object doWork() {
            return null;
        }
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testRsIDFeatureContext(final VariantContext toAnnotate, final List<VariantContext> dbSNPRecords, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureContext ctx = new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(file)), new SimpleInterval(toAnnotate));
        final VariantOverlapAnnotator annotator = makeAnnotator(file, "dbsnp");
        final VariantContext annotated = annotator.annotateRsID(ctx, toAnnotate);
        Assert.assertNotNull(annotated);
        Assert.assertEquals(annotated, toAnnotate);   //nothing at given position
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testRsIDFeatureContextWith2Sets(final VariantContext toAnnotate, final List<VariantContext> dbSNPRecords, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureContext ctx = new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(file)), new SimpleInterval(toAnnotate));
        final VariantOverlapAnnotator annotator = makeAnnotator(file, "dbsnp", "binding");
        final VariantContext annotated = annotator.annotateRsID(ctx, toAnnotate);
        Assert.assertNotNull(annotated);
        Assert.assertEquals(annotated, toAnnotate);   //nothing at given position
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testRsIDFeatureContextWithNoDBSnP(final VariantContext toAnnotate, final List<VariantContext> dbSNPRecords, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureContext ctx = new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(file)), new SimpleInterval(toAnnotate));
        final VariantOverlapAnnotator annotator = makeAnnotator(file, null, "binding");
        final VariantContext annotated = annotator.annotateRsID(ctx, toAnnotate);
        Assert.assertNotNull(annotated);
        Assert.assertEquals(annotated, toAnnotate);   //nothing at given position
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testAnnotateOverlapsFeatureContext(final VariantContext toAnnotate, final List<VariantContext> records, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureContext ctx = new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(file)), new SimpleInterval(toAnnotate));
        final VariantOverlapAnnotator annotator = makeAnnotator(file, "dbsnp");
        final VariantContext annotated = annotator.annotateOverlaps(ctx, toAnnotate);
        Assert.assertEquals(annotated, toAnnotate);   //nothing at given position
    }

    @Test(dataProvider = "AnnotateRsIDData")
    public void testAnnotateOverlapsFeatureContextWith2Sets(final VariantContext toAnnotate, final List<VariantContext> records, final String expectedID, final boolean expectOverlap) throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureContext ctx = new FeatureContext(new FeatureManager(new ArtificialFeatureContainingCommandLineProgram_ForVariantOverlap(file)), new SimpleInterval(toAnnotate));
        final VariantOverlapAnnotator annotator = makeAnnotator(file, "dbsnp", "binding");
        final VariantContext annotated = annotator.annotateOverlaps(ctx, toAnnotate);
        Assert.assertEquals(annotated, toAnnotate);   //nothing at given position
    }
}
