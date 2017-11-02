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
public class WindowStepArgumentCollectionTest extends BaseTest{

    private static class FakeToolWStep {
        @ArgumentCollection
        final WindowStepArgumentCollection ws;

        FakeToolWStep(final WindowStepArgumentCollection ws) {
            this.ws = ws;
        }
    }


    @DataProvider(name="illegalWStep")
    public Object[][] getIllegalSteps() {
        return new Object[][] {{-10}, {-1}, {0}};
    }

    @DataProvider(name="legalWStep")
    public Object[][] getLegalSteps() {
        return new Object[][] {{1}, {5}, {10}};
    }

    @Test(dataProvider = "illegalWStep", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalDefaultWindowSize(final int ws) {
        final WindowStepArgumentCollection ac = new DefaultWindowStepArgumentCollection(ws);
        System.err.println("Window step: "+ac.getWindowStep());
    }

    @Test(dataProvider = "illegalWStep", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalRequiredWindowSize(final int ws) {
        final WindowStepArgumentCollection ac = new OptionalWindowStepArgumentCollection(ws);
        System.err.println("Window step: "+ac.getWindowStep());
    }

    @Test(dataProvider = "legalWStep")
    public void testLegalDefaultWindowSize(final int ws) {
        final WindowStepArgumentCollection ac = new DefaultWindowStepArgumentCollection(ws);
        Assert.assertEquals(ac.getWindowStep(), ws);
    }

    @Test(dataProvider = "legalWStep")
    public void testLegalRequiredWindowSize(final int ws) {
        final WindowStepArgumentCollection ac = new OptionalWindowStepArgumentCollection(ws);
        Assert.assertEquals(ac.getWindowStep(), ws);
    }

    @DataProvider(name="defaultWStep")
    private Object[][] withDefault() {
        return new Object[][]{
                {new FakeToolWStep(new DefaultWindowStepArgumentCollection(10)), 10},
                {new FakeToolWStep(new DefaultWindowStepArgumentCollection(5)),  5}};
    }

    @DataProvider(name="requiredWStep")
    private Object[][] withRequired() {
        return new Object[][]{
                {new FakeToolWStep(new OptionalWindowStepArgumentCollection(10)), 10},
                {new FakeToolWStep(new OptionalWindowStepArgumentCollection(5)), 5}};
    }

    @DataProvider(name="bothWStep")
    private Object[][] getBoth() {
        return ArrayUtils.addAll(withDefault(), withRequired());
    }

    @Test(dataProvider = "bothWStep")
    public void testNoOptionProvided(final FakeToolWStep tool, final int defaultValue){
        String[] args = {};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.ws.getWindowStep(), defaultValue);
    }

    @Test(dataProvider = "requiredWStep")
    public void provideOptionToRequired(final FakeToolWStep tool, final int defaultValue) {
        String[] args = {"--windowStep", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        Assert.assertEquals(tool.ws.getWindowStep(), defaultValue+1);
    }

    @Test(dataProvider = "defaultWStep", expectedExceptions = UserException.CommandLineException.class)
    public void provideOptionToDefault(final FakeToolWStep tool, final int defaultValue) {
        String[] args = {"--windowStep", String.valueOf(defaultValue+1)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
    }

    @Test(dataProvider = "illegalWStep", expectedExceptions = UserException.BadArgumentValue.class)
    public void provideIllegalValue(final int ws) {
        final FakeToolWStep tool = new FakeToolWStep(new OptionalWindowStepArgumentCollection(1000));
        String[] args = {"--windowStep", String.valueOf(ws)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userStep = tool.ws.getWindowStep();
        System.err.println("User window-step: "+userStep);
    }

    @Test(dataProvider = "legalWStep")
    public void provideLegalValue(final int ws) {
        final FakeToolWStep tool = new FakeToolWStep(new OptionalWindowStepArgumentCollection(1000));
        String[] args = {"--windowStep", String.valueOf(ws)};
        CommandLineParser clp = new CommandLineParser(tool);
        clp.parseArguments(System.out, args);
        int userStep = tool.ws.getWindowStep();
        Assert.assertEquals(userStep, ws);
    }
}