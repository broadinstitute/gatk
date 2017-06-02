package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Unit tests for {@link TargetArgumentCollection}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetArgumentCollectionUnitTest {

    private static final File BAD_TEST_TARGETS = new File("/not-a-file.tab");
    private static final File TEST_TARGETS = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/targetargumentcollection-test-targets.tab");

    @Test
    public void testOptionalWithoutDefaultSupplierAndMissingTargetsFile() {
        final TargetArgumentCollection subject = new TargetArgumentCollection();
        Assert.assertNull(subject.readTargetCollection(true));
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testRequiredWithoutDefaultSupplierAndMissingTargetsFile() {
        final TargetArgumentCollection subject = new TargetArgumentCollection();
        subject.readTargetCollection(false);
    }

    @Test
    public void testOptionalWithDefaultNullSupplierAndMissingTargetsFile() {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> null);
        Assert.assertNull(subject.readTargetCollection(true));
    }

    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
    public void testRequiredWithDefaultNullSupplierAndMissingTargetsFile() {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> null);
        subject.readTargetCollection(false);
    }

    @Test
    public void testOptionalWithoutDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection();
        subject.targetsFile = TEST_TARGETS;
        final TargetCollection<Target> targets = subject.readTargetCollection(true);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test
    public void testRequiredWithoutDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection();
        subject.targetsFile = TEST_TARGETS;
        final TargetCollection<Target> targets = subject.readTargetCollection(true);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test
    public void testOptionalWithNullDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> null);
        subject.targetsFile = TEST_TARGETS;
        final TargetCollection<Target> targets = subject.readTargetCollection(false);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test
    public void testRequiredWithNullDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> null);
        subject.targetsFile = TEST_TARGETS;
        final TargetCollection<Target> targets = subject.readTargetCollection(false);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test
    public void testRequiredWithDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> BAD_TEST_TARGETS);
        subject.targetsFile = TEST_TARGETS;
        final TargetCollection<Target> targets = subject.readTargetCollection(false);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test
    public void testNoExplicitRequiredWithDefaultSupplier() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection(() -> TEST_TARGETS);
        final TargetCollection<Target> targets = subject.readTargetCollection(false);
        Assert.assertNotNull(targets);
        final List<Target> expectedTargets = ReadCountCollectionUtils.parse(TEST_TARGETS).targets();
        final List<Target> actualTargets = targets.targets();
        Assert.assertEquals(actualTargets, expectedTargets);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testBadTargetFile() throws IOException {
        final TargetArgumentCollection subject = new TargetArgumentCollection();
        subject.targetsFile = BAD_TEST_TARGETS;
        subject.readTargetCollection(false);
    }
}
