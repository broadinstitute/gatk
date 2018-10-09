package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedHashSet;

public class LocatableFuncotationCreatorUnitTest extends GATKBaseTest {

    @Test
    public void testCreate() {
        final SimpleInterval locatable = new SimpleInterval("1", 1000, 2000);
        final Funcotation guess = LocatableFuncotationCreator.create(locatable,
                Allele.create("C"), "TEST_NAME");
        Assert.assertEquals(guess.getFieldNames(), new LinkedHashSet<>(
                Arrays.asList(LocatableFuncotationCreator.CONTIG_FIELD_NAME,
                        LocatableFuncotationCreator.START_FIELD_NAME,
                        LocatableFuncotationCreator.END_FIELD_NAME
                )));
        Assert.assertEquals(guess.getField(LocatableFuncotationCreator.CONTIG_FIELD_NAME), locatable.getContig());
        Assert.assertEquals(guess.getField(LocatableFuncotationCreator.START_FIELD_NAME), String.valueOf(locatable.getStart()));
        Assert.assertEquals(guess.getField(LocatableFuncotationCreator.END_FIELD_NAME), String.valueOf(locatable.getEnd()));

        Assert.assertEquals(guess.getMetadata(), LocatableFuncotationCreator.METADATA);

    }
}
