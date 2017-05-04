package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.testng.Assert;
import org.testng.annotations.Test;

public class PSTaxInfoTest {

    @Test
    public void testEqualsAndHashcode() {
        final PSTaxInfo infoA = new PSTaxInfo();
        infoA.parent_tax = "1";
        infoA.length = 1000;
        infoA.name = "tax A";
        infoA.rank = "species";
        infoA.ref_names.add("ref|NC_A1");
        infoA.ref_names.add("ref|NC_A2");

        final PSTaxInfo infoB = new PSTaxInfo();
        infoB.parent_tax = "2";
        infoB.length = 1000;
        infoB.name = "tax B";
        infoB.rank = "species";
        infoB.ref_names.add("ref|NC_B1");

        Assert.assertNotEquals(infoA, infoB);
        Assert.assertNotEquals(infoA.hashCode(), infoB.hashCode());

        final PSTaxInfo infoA2 = new PSTaxInfo();
        infoA2.parent_tax = "1";
        infoA2.length = 1000;
        infoA2.name = "tax A";
        infoA2.rank = "species";
        infoA2.ref_names.add("ref|NC_A2");
        infoA2.ref_names.add("ref|NC_A1");

        Assert.assertEquals(infoA, infoA2);
        Assert.assertNotEquals(infoA.hashCode(), infoA2.hashCode());

        final PSTaxInfo infoA3 = new PSTaxInfo();
        infoA3.parent_tax = "1";
        infoA3.length = 1000;
        infoA3.name = "tax A";
        infoA3.rank = "species";
        infoA3.ref_names.add("ref|NC_A1");
        infoA3.ref_names.add("ref|NC_A2");

        Assert.assertEquals(infoA.hashCode(), infoA3.hashCode());
    }

    @Test
    public void testToString() {
        final PSTaxInfo infoA = new PSTaxInfo();
        infoA.parent_tax = "1";
        infoA.length = 1000;
        infoA.name = "tax A";
        infoA.rank = "species";
        infoA.ref_names.add("ref|NC_A1");
        infoA.ref_names.add("ref|NC_A2");

        final String strA = infoA.toString();
        Assert.assertNotNull(strA);
        Assert.assertTrue(strA.contains(infoA.name));
    }
}