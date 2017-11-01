package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Map;

public final class TandemRepeatUnitTest extends GATKBaseTest {

    @Test
    public void testUsingVC() {

        // - [ref] / ATC from 20-20
        final String insLoc = "chr1";
        final int insLocStart = 2;
        final int insLocStop = 2;

        final byte[] refBytes = "GTATCATCATCGGA".getBytes();

        final Allele nullR = Allele.create("A", true);
        final Allele atc   = Allele.create("AATC", false);

        // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
        final VariantContext vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStop, Arrays.asList(nullR,atc)).make();

        // we test that the interval from which the ReferenceContext is constructed does not need to exactly overlap
        // the VariantContext.  The annotation should be able to handle this.
        final SimpleInterval interval= new SimpleInterval(insLoc, insLocStart + 3, insLocStop + 4);

        final SimpleInterval interval1 = new SimpleInterval(insLoc, 1, refBytes.length);
        final ReferenceBases ref1 = new ReferenceBases(refBytes, interval1);

        final SAMSequenceDictionary dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(insLoc, refBytes.length)));
        final ReferenceContext ref = new ReferenceContext(ReferenceDataSource.of(ref1, dict), interval, 20, 20);
        final InfoFieldAnnotation ann = new TandemRepeat();
        final Map<String, Object> a = ann.annotate(ref, vc, null);

        Assert.assertEquals(a.size(), 3);
        Assert.assertEquals(a.get(GATKVCFConstants.STR_PRESENT_KEY), true);
        Assert.assertEquals(a.get(GATKVCFConstants.REPEAT_UNIT_KEY), "ATC");
        Assert.assertEquals(a.get(GATKVCFConstants.REPEATS_PER_ALLELE_KEY), Arrays.asList(3,4)); //3 repeats in the reference, 4 in the variant
    }

    @Test
    public void testUsingVCNoRepeat() {

        // - [ref] / ATC from 20-20
        final String insLoc = "chr1";
        final int insLocStart = 6;
        final int insLocStop = 6;

        final byte[] refBytes = "GTATCATCATCGGA".getBytes();

        final Allele nullR = Allele.create("A", true);
        final Allele atc   = Allele.create("AATC", false);

        // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
        final VariantContext vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStop, Arrays.asList(nullR,atc)).make();

        final SimpleInterval interval= new SimpleInterval(insLoc, insLocStart, insLocStop);

        final SimpleInterval interval1 = new SimpleInterval(insLoc, 1, refBytes.length);
        final ReferenceBases ref1 = new ReferenceBases(refBytes, interval1);

        final SAMSequenceDictionary dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(insLoc, refBytes.length)));
        final ReferenceContext ref = new ReferenceContext(ReferenceDataSource.of(ref1, dict), interval, 0, 20);

        final InfoFieldAnnotation ann = new TandemRepeat();
        final Map<String, Object> a = ann.annotate(ref, vc, null);

        Assert.assertTrue(a.isEmpty());
    }

    @Test
    public void testUsingVCNotIndel() {

        // - [ref] / ATC from 20-20
        String insLoc = "chr1";
        int insLocStart = 2;
        int insLocStop = 2;

        byte[] refBytes = "GTATCATCATCGGA".getBytes();

        Allele nullR = Allele.create("A", true);
        Allele atc   = Allele.create("C", false);

        // A*,ATC, context = ATC ATC ATC : (ATC)3 -> (ATC)4
        VariantContext vc = new VariantContextBuilder("foo", insLoc, insLocStart, insLocStop, Arrays.asList(nullR,atc)).make();

        final SimpleInterval interval= new SimpleInterval("chr1", insLocStart, insLocStop);
        final String contigName = "chr1";
        final SimpleInterval interval1 = new SimpleInterval(contigName, 1, refBytes.length);
        final ReferenceBases ref1 = new ReferenceBases(refBytes, interval1);

        final SAMSequenceDictionary dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(contigName, refBytes.length)));
        final ReferenceContext ref = new ReferenceContext(ReferenceDataSource.of(ref1, dict), interval, 0, 20);
        final InfoFieldAnnotation ann = new TandemRepeat();
        final Map<String, Object> a = ann.annotate(ref, vc, null);

        Assert.assertTrue(a.isEmpty());
    }

    @Test
    public void testDescriptions() {
        final InfoFieldAnnotation ann = new TandemRepeat();
        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STR_PRESENT_KEY),
                GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.REPEAT_UNIT_KEY),
                GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.REPEATS_PER_ALLELE_KEY)
        ));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(
                GATKVCFConstants.STR_PRESENT_KEY,
                GATKVCFConstants.REPEAT_UNIT_KEY,
                GATKVCFConstants.REPEATS_PER_ALLELE_KEY));
    }
}
