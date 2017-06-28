package org.broadinstitute.hellbender.engine.filters;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class CigarReadFilterTest {

    @DataProvider
    public Object[][] createCigarFilterStrings(){

        return new Object[][] {
                {"^", false},
                {"$", false},
                {"*", true },
                {"H", true},
                {"HS", true},
                {"HHHHSSSSDSSSSHHH", true},
                {"SH", true},
                {"SHS", false},
                {"SHSH", false},
                {"=H", true},
                {"105H", true},
                {"329H3S", true},
                {"329", false},
                {"^M$", true},
                {"M", true},
                {"M%", false},
                {"M*", true},
                {"*M", true},
                {"10P", true},
                {"=10P", true},
                {"<=P", false},
                {"<=9=", true},
                {"10P<>34X", false},
                {"10P<>=34X", false},
                {"10P=34X", true},
                {"10P<34X", true},
                {"10P<=34X", true},
                {"10P>34X", true},
                {"10P>=34X", true},
                {"10P>=34X3S2H", true},
                {"12H<3S10P>=34X3S2H", true},
                {"12H<3SH10P>=34X3S2H", false},
                {"^12H<3S10P>=34X3S2H$", true},
        };
    }

    @Test(dataProvider="createCigarFilterStrings")
    public void testValidate(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.validate(cigarFilterString), isValid );
    }

    @DataProvider
    public Object[][] createCigarMatchElementStrings() {

//        "(\\^?(?:[<>]|[<>]=)?\\d*[SHMIDNPX=*]\\$?)"
        return new Object[][] {
                {"^", false},
                {"$", false},
                {"M$", true},
                {"^*$", true},
                {"107X$", true},
                {"=19X", false},
                {"<=19X", true},
                {">=19X", true},
                {">19X", true},
                {"<19X", true},
        };
    }

    @Test(dataProvider="createCigarMatchElementStrings")
    public void testNextCigarMatchElementPattern(String cigarFilterString, boolean isValid) {
        Assert.assertEquals( CigarReadFilter.nextCigarMatchElementPattern.matcher(cigarFilterString).matches(), isValid );
    }

//    @Test
//    public void testSetDescription() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
//
//    @Test
//    public void testReadFilterTestMethod() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
//
//    @Test
//    public void testConstructor() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
//
//    @Test
//    public void testGetDescription() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
//
//    @Test
//    public void testGetMatchElementCollection() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
//
//    @Test
//    public void testCreateMatchElementFromMatchString() {
//        //FIXME
//        throw new RuntimeException("MUST BE IMPLEMENTED");
//    }
}
