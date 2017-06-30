package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.ImmutableMap;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link GermlinePloidyAnnotatedTargetCollection}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class GermlinePloidyAnnotatedTargetCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = TestResources.publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File CONTIG_PLOIDY_ANNOTS_FILE = new File(TEST_SUB_DIR, "contig_annots.tsv");
    private static final File TARGET_FILE = new File(TEST_SUB_DIR, "sex_genotyper_agilent_targets_trunc.tsv");

    private static List<ContigGermlinePloidyAnnotation> contigPloidyAnnotsList;
    private static List<Target> targetList;
    private static GermlinePloidyAnnotatedTargetCollection ploidyAnnotsCollection;

    @BeforeClass
    public static void loadTestResources() {
        contigPloidyAnnotsList = ContigGermlinePloidyAnnotationTableReader.readContigGermlinePloidyAnnotationsFromFile(CONTIG_PLOIDY_ANNOTS_FILE);
        targetList = TargetTableReader.readTargetFile(TARGET_FILE);
        ploidyAnnotsCollection = new GermlinePloidyAnnotatedTargetCollection(contigPloidyAnnotsList, targetList);
    }

    @Test
    public void testTargetPloidy() {
        /* on an 1 contig */
        final Target testAutosomalTarget = new Target("UNNAMED_TARGET_0");
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testAutosomalTarget, "SEX_XX"), 2);
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testAutosomalTarget, "SEX_XX"), 2);

        /* on X contig */
        final Target testXTarget = new Target("UNNAMED_TARGET_10");
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testXTarget, "SEX_XX"), 2);
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testXTarget, "SEX_XY"), 1);

        /* on Y contig */
        final Target testYTarget = new Target("UNNAMED_TARGET_19");
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testYTarget, "SEX_XX"), 0);
        Assert.assertEquals(ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(testYTarget, "SEX_XY"), 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNonexistentPloidyTag() {
        ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(new Target("UNNAMED_TARGET_0"), "FOO");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNonexistentTarget() {
        ploidyAnnotsCollection.getTargetGermlinePloidyByGenotype(new Target("UNNAMED_TARGET_FOO"), "SEX_XX");
    }

    @Test
    public void testAutosomalTargetsList() {
        final Set<String> autosomalTargetNamesResult = ploidyAnnotsCollection.getAutosomalTargetList()
                .stream().map(Target::getName).collect(Collectors.toSet());
        final Set<String> autosomalTargetNamesExpected = Arrays.stream(new String[] {
                "UNNAMED_TARGET_0", "UNNAMED_TARGET_1", "UNNAMED_TARGET_2",
                "UNNAMED_TARGET_3", "UNNAMED_TARGET_4", "UNNAMED_TARGET_5",
                "UNNAMED_TARGET_6", "UNNAMED_TARGET_7", "UNNAMED_TARGET_8",
                "UNNAMED_TARGET_9"
        }).collect(Collectors.toSet());
        Assert.assertTrue(autosomalTargetNamesResult.equals(autosomalTargetNamesExpected));
    }

    @Test
    public void testAllosomalTargetsList() {
        final Set<String> allosomalTargetNamesResult = ploidyAnnotsCollection.getAllosomalTargetList()
                .stream().map(Target::getName).collect(Collectors.toSet());
        final Set<String> allosomalTargetNamesExpected = Arrays.stream(new String[] {
                "UNNAMED_TARGET_10", "UNNAMED_TARGET_11", "UNNAMED_TARGET_12",
                "UNNAMED_TARGET_13", "UNNAMED_TARGET_14", "UNNAMED_TARGET_15",
                "UNNAMED_TARGET_16", "UNNAMED_TARGET_17", "UNNAMED_TARGET_18",
                "UNNAMED_TARGET_19", "UNNAMED_TARGET_20", "UNNAMED_TARGET_21"
        }).collect(Collectors.toSet());
        Assert.assertTrue(allosomalTargetNamesResult.equals(allosomalTargetNamesExpected));
    }

    /* check if targets have unique names */
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNonUniqueTargetNames() {
        final List<Target> badTargetList = Arrays.stream(new Target[] {
                new Target("NON_UNIQUE_NAME", new SimpleInterval("1", 1, 10)),
                new Target("NON_UNIQUE_NAME", new SimpleInterval("1", 10, 20))
        }).collect(Collectors.toList());
        new GermlinePloidyAnnotatedTargetCollection(contigPloidyAnnotsList, badTargetList);
    }

    /* check if ploidy map is empty */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testEmptyPloidyMap() {
        final List<ContigGermlinePloidyAnnotation> badContigAnnotsList = Arrays.stream(new ContigGermlinePloidyAnnotation[] {
                new ContigGermlinePloidyAnnotation("1", ContigClass.AUTOSOMAL, new HashMap<>())
        }).collect(Collectors.toList());
        new GermlinePloidyAnnotatedTargetCollection(badContigAnnotsList, targetList);
    }
    /* check if annotated contigs are unique */
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMultiplyAnnotatedContigs() {
        final Map<String, Integer> ploidyMap = ImmutableMap.<String, Integer>builder()
                .put("SEX_XX", 2).put("SEX_XY", 2).build();
        final List<ContigGermlinePloidyAnnotation> badContigAnnotsList = Arrays.stream(new ContigGermlinePloidyAnnotation[] {
                new ContigGermlinePloidyAnnotation("1", ContigClass.AUTOSOMAL, ploidyMap),
                new ContigGermlinePloidyAnnotation("1", ContigClass.ALLOSOMAL, ploidyMap)
        }).collect(Collectors.toList());
        new GermlinePloidyAnnotatedTargetCollection(badContigAnnotsList, targetList);
    }

    /* check if all contigs are annotated */
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingAnnotations() {
        final Map<String, Integer> ploidyMap = ImmutableMap.<String, Integer>builder()
                .put("SEX_XX", 2).put("SEX_XY", 2).build();
        final List<ContigGermlinePloidyAnnotation> badContigAnnotsList = Arrays.stream(new ContigGermlinePloidyAnnotation[] {
                new ContigGermlinePloidyAnnotation("1", ContigClass.AUTOSOMAL, ploidyMap)
        }).collect(Collectors.toList());
        new GermlinePloidyAnnotatedTargetCollection(badContigAnnotsList, targetList);
    }
}
