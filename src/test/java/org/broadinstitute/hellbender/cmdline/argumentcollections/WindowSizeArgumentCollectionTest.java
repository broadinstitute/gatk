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
public class WindowSizeArgumentCollectionTest extends BaseTest{

    private static class FakeToolWSize {
        @ArgumentCollection
        final WindowSizeArgumentCollection ws;

        FakeToolWSize(final WindowSizeArgumentCollection ws) {
            this.ws = ws;
        }
    }

    @DataProvider(name="illegalWSize")
    public Object[][] getIllegalSizes() {
        return new Object[][] {{-10}, {-1}, {0}};
    }

    @DataProvider(name="legalWSize")
    public Object[][] getLegalSizes() {
        return new Object[][] {{1}, {5}, {10}};
    }

    @Test(dataProvider = "illegalWSize", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalDefaultWindowSize(final int ws) {
        final WindowSizeArgumentCollection ac = new DefaultWindowSizeArgumentCollection(ws);
        System.err.println("Window size: "+ac.getWindowSize());
    }

    @Test(dataProvider = "illegalWSize", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalRequiredWindowSize(final int ws) {
        final WindowSizeArgumentCollection ac = new OptionalWindowSizeArgumentCollection(ws);
        System.err.println("Window size: "+ac.getWindowSize());
    }

    @Test(dataProvider = "legalWSize")
    public void testLegalDefaultWindowSize(final int ws) {
        final WindowSizeArgumentCollection ac = new DefaultWindowSizeArgumentCollection(ws);
        Assert.assertEquals(ac.getWindowSize(), ws);
    }

    @Test(dataProvider = "legalWSize")
    public void testLegalRequiredWindowSize(final int ws) {
        final WindowSizeArgumentCollection ac = new OptionalWindowSizeArgumentCollection(ws);
        Assert.assertEquals(ac.getWindowSize(), ws);
    }

    @DataProvider(name="defaultWPadding")
    private Object[][] withDefault() {
        return new Object[][]{
                {new FakeToolWSize(new DefaultWindowSizeArgumentCollection(10)), 10},
                {new FakeToolWSize(new DefaultWindowSizeArgumentCollection(5)),  5}};
    }

    @DataProvider(name="requiredWPadding")
    private Object[][] withRequired() {
        return new Object[][]{
                {new FakeToolWSize(new OptionalWindowSizeArgumentCollection(10)), 10},
                {new FakeToolWSize(new OptionalWindowSizeArgumentCollection(5)),   5}};
    }

    @DataProvider(name="bothWPadding")
    private Object[][] getBoth() {
        return ArrayUtils.addAll(withDefault(), withRequired());
    }

    @Test(dataProvider = "bothWPadding")
    public void testNoOptionProvided(final FakeToolWSize tool, final int defaultValue){
        String[] args = {};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.ws.getWindowSize(), defaultValue);
    }

    @Test(dataProvider = "requiredWPadding")
    public void provideOptionToRequired(final FakeToolWSize tool, final int defaultValue) {
        String[] args = {"--windowSize", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.ws.getWindowSize(), defaultValue+1);
    }

    @Test(dataProvider = "defaultWPadding", expectedExceptions = UserException.CommandLineException.class)
    public void provideOptionToDefault(final FakeToolWSize tool, final int defaultValue) {
        String[] args = {"--windowSize", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
    }

    @Test(dataProvider = "illegalWSize", expectedExceptions = UserException.BadArgumentValue.class)
    public void provideIllegalValue(final int ws) {
        final FakeToolWSize tool = new FakeToolWSize(new OptionalWindowSizeArgumentCollection(1000));
        String[] args = {"--windowSize", String.valueOf(ws)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userSize = tool.ws.getWindowSize();
        System.err.println("User window-size: "+userSize);
    }

    @Test(dataProvider = "legalWSize")
    public void provideLegalValue(final int ws) {
        final FakeToolWSize tool = new FakeToolWSize(new OptionalWindowSizeArgumentCollection(1000));
        String[] args = {"--windowSize", String.valueOf(ws)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userSize = tool.ws.getWindowSize();
        Assert.assertEquals(userSize, ws);
    }
}