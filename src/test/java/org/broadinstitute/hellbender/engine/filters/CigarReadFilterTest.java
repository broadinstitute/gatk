package org.broadinstitute.hellbender.engine.filters;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class CigarReadFilterTest {

    @DataProvider
    public Object[][] createCigarFilterStrings(){

        //"\\^*$|\\^(?:\\d*H)?(?:\\d*S)?(?:\\d*[MIDNPX=%])*(?:\\d*S)*(?:\\d*H)*\\$"
        return new Object[][]{
                {"^", false},
                {"$", false},
                {"*", true },
                {"H", true},
                {"HS", true},
                {"SH", true},
                {"SHSH", false},
                {"^", false},
                {"^", false},
                {"^", false},
                {"^", false},
        };
    }

    @Test(dataProvider="createCigarFilterStrings")
    public void testValidate(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), isValid );
    }

}
