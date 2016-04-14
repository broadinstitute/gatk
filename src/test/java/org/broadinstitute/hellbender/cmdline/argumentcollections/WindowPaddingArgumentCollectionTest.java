package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineParser;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class WindowPaddingArgumentCollectionTest extends BaseTest{

    private static class FakeToolWPadding {
        @ArgumentCollection
        final WindowPaddingArgumentCollection wp;

        FakeToolWPadding(final WindowPaddingArgumentCollection wp) {
            this.wp = wp;
        }
    }

    @DataProvider(name="illegalWPadding")
    public Object[][] getIllegalPadding() {
        return new Object[][] {{-10}, {-5}, {-1}};
    }

    @DataProvider(name="legalWPadding")
    public Object[][] getLegalPadding() {
        return new Object[][] {{0}, {1}, {5}};
    }

    @Test(dataProvider = "illegalWPadding", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalDefaultWindowSize(final int pad) {
        final WindowPaddingArgumentCollection ac = new DefaultWindowPaddingArgumentCollection(pad);
        System.err.println("Window padding: "+ac.getWindowPadding());
    }

    @Test(dataProvider = "illegalWPadding", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalRequiredWindowSize(final int pad) {
        final WindowPaddingArgumentCollection ac = new OptionalWindowPaddingArgumentCollection(pad);
        System.err.println("Window padding: "+ac.getWindowPadding());
    }

    @Test(dataProvider = "legalWPadding")
    public void testLegalDefaultWindowSize(final int pad) {
        final WindowPaddingArgumentCollection ac = new DefaultWindowPaddingArgumentCollection(pad);
        Assert.assertEquals(ac.getWindowPadding(), pad);
    }

    @Test(dataProvider = "legalWPadding")
    public void testLegalRequiredWindowSize(final int pad) {
        final WindowPaddingArgumentCollection ac = new OptionalWindowPaddingArgumentCollection(pad);
        Assert.assertEquals(ac.getWindowPadding(), pad);
    }

    @DataProvider(name="defaultWPadding")
    private Object[][] withDefault() {
        return new Object[][]{
                {new FakeToolWPadding(new DefaultWindowPaddingArgumentCollection(10)), 10},
                {new FakeToolWPadding(new DefaultWindowPaddingArgumentCollection(5)),  5}};
    }

    @DataProvider(name="requiredWPadding")
    private Object[][] withRequired() {
        return new Object[][]{
                {new FakeToolWPadding(new OptionalWindowPaddingArgumentCollection(10)), 10},
                {new FakeToolWPadding(new OptionalWindowPaddingArgumentCollection(5)),   5}};
    }

    @DataProvider(name="bothWPadding")
    private Object[][] getBoth() {
        return ArrayUtils.addAll(withDefault(), withRequired());
    }

    @Test(dataProvider = "bothWPadding")
    public void testNoOptionProvided(final FakeToolWPadding tool, final int defaultValue){
        String[] args = {};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.wp.getWindowPadding(), defaultValue);
    }

    @Test(dataProvider = "requiredWPadding")
    public void provideOptionToRequired(final FakeToolWPadding tool, final int defaultValue) {
        String[] args = {"--windowPadding", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.wp.getWindowPadding(), defaultValue+1);
    }

    @Test(dataProvider = "defaultWPadding", expectedExceptions = UserException.CommandLineException.class)
    public void provideOptionToDefault(final FakeToolWPadding tool, final int defaultValue) {
        String[] args = {"--windowPadding", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
    }

    @Test(dataProvider = "illegalWPadding", expectedExceptions = UserException.BadArgumentValue.class)
    public void provideIllegalValue(final int pad) {
        final FakeToolWPadding tool = new FakeToolWPadding(new OptionalWindowPaddingArgumentCollection(1000));
        String[] args = {"--windowPadding", String.valueOf(pad)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userPadding = tool.wp.getWindowPadding();
        System.err.println("User padding: "+userPadding);
    }

    @Test(dataProvider = "legalWPadding")
    public void provideLegalValue(final int pad) {
        final FakeToolWPadding tool = new FakeToolWPadding(new OptionalWindowPaddingArgumentCollection(1000));
        String[] args = {"--windowPadding", String.valueOf(pad)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userPadding = tool.wp.getWindowPadding();
        Assert.assertEquals(userPadding, pad);
    }
}