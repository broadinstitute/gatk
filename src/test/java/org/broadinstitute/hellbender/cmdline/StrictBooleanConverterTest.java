package org.broadinstitute.hellbender.cmdline;

import joptsimple.ValueConversionException;
import org.testng.Assert;
import org.testng.annotations.Test;


public final class StrictBooleanConverterTest {
    @Test
    public void recognizedValues(){
        StrictBooleanConverter converter = new StrictBooleanConverter();
        Assert.assertEquals("true", converter.convert("true"));
        Assert.assertEquals("true", converter.convert("T"));
        Assert.assertEquals("true",converter.convert("TRUE"));
        Assert.assertEquals("false",converter.convert("F"));
        Assert.assertEquals("false", converter.convert("False"));
    }

    @Test(expectedExceptions = ValueConversionException.class)
    public void unrecognizedValues(){
        StrictBooleanConverter converter = new StrictBooleanConverter();
        converter.convert("unprovable");
    }

}
