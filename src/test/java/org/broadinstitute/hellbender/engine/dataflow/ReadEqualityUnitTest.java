package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Position;
import com.google.api.services.genomics.model.LinearAlignment;
import com.google.api.services.genomics.model.Read;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;


public class ReadEqualityUnitTest extends BaseTest {

    // Test to prove that Google Genomics Reads have a working equals() method
    @Test
    public void testEquality() {
        Read read1 = new Read();
        Read read2 = new Read();

        read1.setAlignment(new LinearAlignment());
        read1.getAlignment().setPosition(new Position());
        read1.getAlignment().getPosition().setReferenceName("FOO");

        read2.setAlignment(new LinearAlignment());
        read2.getAlignment().setPosition(new Position());
        read2.getAlignment().getPosition().setReferenceName("FOO");

        Assert.assertEquals(read1, read2, "equal reads not equal");

        read2.getAlignment().getPosition().setReferenceName("BAR");
        Assert.assertNotEquals(read1, read2, "unequal reads are equal");
    }
}
