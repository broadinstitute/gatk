package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.io.output.NullOutputStream;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.mockito.internal.util.collections.Sets;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.collections.Lists;

import java.io.PrintStream;
import java.util.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GATKAnnotationPluginDescriptorUnitTest extends GATKBaseTest {
    // null print stream for the tests
    private static final PrintStream nullMessageStream = new PrintStream(new NullOutputStream());

    private GATKAnnotationArgumentCollection getEmptyArgs() {
        return new DefaultGATKVariantAnnotationArgumentCollection();
    }

//======================================================================================================================
// Methods for computing individual annotations
    private static final Allele Aref = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele C = Allele.create("C");

    //make sure that compound hets (with no ref) don't add to het count
    VariantContext inbreedingCoefficientVC = makeVC("1", Arrays.asList(Aref, T, C),
            makeG("s1", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
            makeG("s4", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s6", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s8", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s10", Aref, T, 2530, 0, 7099, 366, 3056, 14931),

            //add a bunch of hom samples that will be ignored if we use s1..s10 as founders
            makeG("s11", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s12", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s13", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s14", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s15", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s16", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931),
            makeG("s17", T, T, 7099, 2530, 0, 7099, 366, 3056, 14931)
    );
    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... genotypes) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(Arrays.asList(genotypes)).filters((Set<String>) null).make();
    }
    private Genotype makeG(String sample, Allele a1, Allele a2, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }
//======================================================================================================================

    private List<Annotation> instantiateAnnotations(final CommandLineParser clp) {
        GATKAnnotationPluginDescriptor annotationPlugin = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        return annotationPlugin.getResolvedInstances();
    }

    @DataProvider(name = "defaultAnnotationsForAllowedValues")
    public Object[][] defaultAnnotationsForAllowedValues() {
        return new Object[][] {
                {Collections.emptyList()},
                {Collections.singletonList(new InbreedingCoeff())},
                {Arrays.asList(new RMSMappingQuality(), new InbreedingCoeff())},
                {Arrays.asList(new InbreedingCoeff(), new RMSMappingQuality())}
        };
    }

    @Test(dataProvider = "defaultAnnotationsForAllowedValues")
    public void testGetAllowedValuesForDescriptorHelp(final List<Annotation> defaultAnnotations) {
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), defaultAnnotations, null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[]{});

        final GATKAnnotationPluginDescriptor pluginDescriptor = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);

        // the help for exclude-annotation should point out to the default Annotations in the order provided
        Assert.assertEquals(pluginDescriptor.getAllowedValuesForDescriptorHelp(StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME),
                defaultAnnotations.stream().map(rf -> rf.getClass().getSimpleName()).collect(Collectors.toSet()));

        // test if the help for annotation is not empty after parsing: if custom validation throws, the help should print the annotations available
        // the complete set could not checked because that requires to discover all the implemented annotations
        Assert.assertFalse(pluginDescriptor.getAllowedValuesForDescriptorHelp(StandardArgumentDefinitions.ANNOTATION_LONG_NAME).isEmpty());

        // test if the help for annotation group is not empty after parsing: if custom validation throws, the help should print the annotations available
        // the complete set could not checked because that requires knowing all the implemented annotation groups
        Assert.assertFalse(pluginDescriptor.getAllowedValuesForDescriptorHelp(StandardArgumentDefinitions.ANNOTATION_GROUP_LONG_NAME).isEmpty());
    }

    @Test
    public void testGetAllowedValuesForDescriptorHelpGroups() {
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, Collections.singletonList(AS_StandardAnnotation.class))),
        Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[]{});

        final GATKAnnotationPluginDescriptor pluginDescriptor = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);

        // the help for exclude-annotation should include annotations that were added to the descriptor by default groups
        Assert.assertTrue(pluginDescriptor.getAllowedValuesForDescriptorHelp(StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME).contains(AS_RMSMappingQuality.class.getSimpleName()));
    }

    @Test
    public void testGetDefaultInstances() {
        final CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(),  Collections.singletonList(RMSMappingQuality.getInstance()), Collections.singletonList(AS_StandardAnnotation.class))),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, new String[]{});

        final GATKAnnotationPluginDescriptor pluginDescriptor = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);

        // Assert that we see the annotation we asked for explicitly
        Assert.assertTrue(pluginDescriptor.getDefaultInstances().stream().anyMatch(a -> a instanceof RMSMappingQuality));
        // Assert but not one we didn't
        Assert.assertFalse(pluginDescriptor.getDefaultInstances().stream().anyMatch(a -> a instanceof ExcessHet));
        // But still see the group annotations
        Assert.assertTrue(pluginDescriptor.getDefaultInstances().stream().anyMatch(a -> a instanceof AS_RMSMappingQuality));;
    }

    @DataProvider
    public Object[][] badAnnotationGroupsDataProvider() {
        Object[][] out = {
                { Arrays.asList(RMSMappingQuality.class)},
                { Arrays.asList(StandardAnnotation.class, Annotation.class)},
                { Arrays.asList(Object.class)}};
        return out;
    }

    @Test (dataProvider = "badAnnotationGroupsDataProvider", expectedExceptions = GATKException.class)
    public void testInvalidToolDefaultAnnotationGroups(List<Class<? extends Annotation>> testGroups) {
        //This test asserts that the plugin descriptor will throw if an invalid annotation group is requested
        new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, testGroups)),
                Collections.emptySet());
    }

    @Test (dataProvider = "badAnnotationGroupsDataProvider", expectedExceptions = CommandLineException.class)
    public void testInvalidRequestedAnnotationGroups(List<Class<? extends Annotation>> testGroups) {
        //This test asserts that the plugin descriptor will throw if an invalid annotation group is requested
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());
        List<String> args = new ArrayList<>();
        for (Class<? extends Annotation> group : testGroups) {
            args.add(StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME);
            args.add(group.getSimpleName());
        }
        clp.parseArguments(nullMessageStream, args.toArray(new String[0]));
    }

    @DataProvider
    public Object[][] badAnnotationArgumentsDataProvider() {
        Object[][] out = {
                { Arrays.asList("-A", "StandardAnnotation")},
                { Arrays.asList("-A", "RMSMappingQuality", "-A", "RMSMappingQuality")},
                { Arrays.asList("-A", "RMSMappingQuality", "-AX", "RMSMappingQuality")},
                { Arrays.asList(StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "StandardAnnotation", "-AX", "RMSMappingQuality", "-AX", "RMSMappingQuality")},
                // { Arrays.asList("-AX", "RMSMappingQuality")}, This is just a warning for now
                { Arrays.asList("-A", "foo")},
                { Arrays.asList("-A", "VariantAnnotator")}};
        return out;
    }

    @Test (dataProvider = "badAnnotationArgumentsDataProvider", expectedExceptions = CommandLineException.class)
    public void testInvalidUserSpecifiedAnnotations(List<String> arguments) {
        //This test asserts that the plugin descriptor will crash if an invalid annotation group is requested
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());
        String[] args = arguments.toArray(new String[arguments.size()]);
        clp.parseArguments(nullMessageStream, args);
    }

    @DataProvider
    public Object[][] annotationsWithArguments(){
        return new Object[][]{{ InbreedingCoeff.class.getSimpleName(), "--founder-id", "s1"}};
    }

    // fail if an annotation with required arguments is specified without corresponding arguments
    // TODO this test is disabled as there are currently no annotations for which the argument is required
    // TODO if any developer does add dependant arguments in the future they should be tested using this format
    @Test(dataProvider = "annotationsWithArguments", expectedExceptions = CommandLineException.MissingArgument.class, enabled = false)
    public void testDependentAnnotationArguments(
            final String annot,
            final String argName,   //unused
            final String[] argValue) { //unused
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());

        String[] args = {"--annotation", annot};  // no args, just enable annotations
        clp.parseArguments(nullMessageStream, args);
    }

    // fail if a annotation's arguments are passed but the annotation itself is not enabled
    @Test(dataProvider = "annotationsWithArguments", expectedExceptions = CommandLineException.class)
    public void testDanglingAnnotationArguments(
            final String annot,
            final String argName,
            final String argValue) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());

        String[] args = { argName, argValue }; // no annotation set

        // no need to instantiate the annotation - dependency errors are caught by the command line parser
        clp.parseArguments(nullMessageStream, args);
    }

    @Test
    public void testAnnotationArguments() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());

        List<String> args = Stream.of("--annotation", InbreedingCoeff.class.getSimpleName() ).collect(Collectors.toList());
        Arrays.asList(new String[]{"s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"}).forEach(arg -> {
            args.addAll(Arrays.asList("--founder-id", arg));
        });

        clp.parseArguments(nullMessageStream, args.toArray(new String[args.size()]));
        List<Annotation> rf = instantiateAnnotations(clp);
        Assert.assertEquals(rf.size(), 1);
        Assert.assertEquals(Double.valueOf((String) ((InbreedingCoeff) rf.get(0))
                        .annotate(null, inbreedingCoefficientVC, null)
                        .get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                -0.3333333, 0.001, "InbreedingCoefficientScores");
    }

    @DataProvider
    public Object[][] artificialAnnotationsTests() {
        return new Object[][]{
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation()),
                        new String[0],
                        true,
                        5,
                        0},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(), new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[0],
                        true,
                        5,
                        -0.3333},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(), new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[]{"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS},
                        false,
                        0,
                        0},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(), new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[]{"--founder-id", "s1","--founder-id", "s2","--founder-id", "s3","--founder-id", "s4","--founder-id", "s5",
                                "--founder-id", "s6","--founder-id", "s7","--founder-id", "s8","--founder-id", "s9","--founder-id", "s10","--founder-id", "s11"},
                        true,
                        5,
                        -0.2941},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(),new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[]{"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-A", "InbreedingCoeff"},
                        false,
                        0,
                        -0.3333},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(),new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[]{"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "StandardAnnotation"},
                        false,
                        0,
                        -0.3333},
                {Arrays.asList(new testChildAnnotation(), new testParentAnnotation(),new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"))),
                        new String[]{"--"+StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "StandardAnnotation", "--founder-id",
                                "s1","--founder-id", "s2","--founder-id", "s3","--founder-id", "s4","--founder-id", "s5", "--founder-id", "s6",
                                "--founder-id", "s7","--founder-id", "s8","--founder-id", "s9","--founder-id", "s10","--founder-id", "s11", "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "ParentAnnotationGroup" },
                        true,
                        5,
                        -0.2941},
                //Testing that pedigree files and founder-ids are additive (s1-s5 are provided by the ped file)
                {Arrays.asList(new InbreedingCoeff()),
                        new String[]{ "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "StandardAnnotation", "-ped",
                                toolsTestDir+"walkers/annotator/testPedigree.ped","--founder-id", "s6",
                                "--founder-id", "s7","--founder-id", "s8","--founder-id", "s9","--founder-id", "s10", "--founder-id", "s11"},
                        false,
                        0,
                        -0.2941},
                //Testing that pedigree files and founder-ids behave when duplicates are supplied (fewer than 10 unique founder-ids between the ped and non-ped so it won't compute)
                {Arrays.asList(new InbreedingCoeff()),
                        new String[]{ "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, "StandardAnnotation", "-ped",
                                toolsTestDir+"walkers/annotator/testPedigree.ped","--founder-id", "s5",
                                "--founder-id", "s4","--founder-id", "s3","--founder-id", "s2","--founder-id", "s1", "--founder-id", "s11"},
                        false,
                        0,
                        0},};


    }

    @Test (dataProvider = "artificialAnnotationsTests")
    public void testMultipleOptionalArguments(final List<Annotation> toolDefaultAnnotations, final String[] arguments, final boolean expectingParent, final int expectedChild, final double inbreedingValue) {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), toolDefaultAnnotations, null)),
                Collections.emptySet());

        clp.parseArguments(nullMessageStream, arguments);
        VariantAnnotatorEngine vae = new VariantAnnotatorEngine(instantiateAnnotations(clp), null, Collections.emptyList(), false);
        VariantContext vc = inbreedingCoefficientVC;
        vc = vae.annotateContext(vc, new FeatureContext(), null, null, a->true);

        Assert.assertEquals(vc.getAttribute("Parent"), expectingParent?"foo":null);
        Assert.assertEquals(vc.getAttributeAsInt("Child", 0), expectedChild);
        Assert.assertEquals(vc.getAttributeAsDouble(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY, 0), inbreedingValue);
    }

    @Test
    public void testHierarchicalAnnotationDiscovery() throws IllegalAccessException, InstantiationException {
        List<Class<? extends Annotation>> toolDefaultAnnotationGroups = Arrays.asList(InbreedingCoeff.class, RMSMappingQuality.class, testChildAnnotation.class, testParentAnnotation.class);
        GATKAnnotationPluginDescriptor pluginDescriptor = new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, Collections.singletonList(ParentAnnotationGroup.class));

        toolDefaultAnnotationGroups.forEach(a -> {
            try {
                pluginDescriptor.createInstanceForPlugin(a);
            } catch (IllegalAccessException | InstantiationException e) {
                e.printStackTrace();
                return;
            }
        });
        pluginDescriptor.validateAndResolvePlugins();

        VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Arrays.asList(pluginDescriptor.getResolvedInstances().toArray(new Annotation[0])), null, Collections.emptyList(), false);
        VariantContext vc = inbreedingCoefficientVC;
        vc = vae.annotateContext(vc, new FeatureContext(), null, null, a->true);

        Assert.assertEquals(vc.getAttribute("Parent"), "foo");
        Assert.assertEquals(vc.getAttributeAsInt("Child", 0), 5);
    }

    @Test
    // Note: this test takes advantage of the fact that InbreedingCoefficient does not compute for pedigrees with fewer than 10 founders.
    public void testToolDefaultAnnotationArgumentsOverridingFromCommandLine() {
        String argName = "--founder-id";
        String[] goodArguments = new String[]{"s1", "s2",  "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"};

        CommandLineParser clp1 = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), Collections.singletonList(new InbreedingCoeff(Collections.emptySet())), null)),
                Collections.emptySet());
        CommandLineParser clp2 = new CommandLineArgumentParser(
                new Object(),
                // Adding a default value which should result in different annotations
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), Collections.singletonList(new InbreedingCoeff(new HashSet<>(Arrays.asList(goodArguments)))), null)),
                Collections.emptySet());

        List<String> completeArguments = Lists.newArrayList("--annotation", InbreedingCoeff.class.getSimpleName());
        List<String> incompleteArguments = Lists.newArrayList("--annotation", InbreedingCoeff.class.getSimpleName(), argName, "s1"); // These arguments are "incomplete" because InbreedingCoefficient
        Arrays.asList(goodArguments).forEach(arg -> {completeArguments.addAll(Arrays.asList(argName, arg));});

        clp1.parseArguments(nullMessageStream, completeArguments.toArray(new String[completeArguments.size()]));
        List<Annotation> goodArgumentsOverridingincompleteArguments = instantiateAnnotations(clp1);
        clp2.parseArguments(nullMessageStream, incompleteArguments.toArray(new String[incompleteArguments.size()]));
        List<Annotation> incompleteArgumentsOverridingGoodArguments = instantiateAnnotations(clp2);

        Assert.assertEquals(goodArgumentsOverridingincompleteArguments.size(), 1);
        Assert.assertEquals(incompleteArgumentsOverridingGoodArguments.size(), 1);
        // assert that with good arguments overriding that the inbreeding coefficient is computed correctly
        Assert.assertEquals(Double.valueOf((String) ((InbreedingCoeff) goodArgumentsOverridingincompleteArguments.get(0))
                        .annotate(null, inbreedingCoefficientVC, null)
                        .get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                -0.3333333, 0.001, "InbreedingCoefficientScores");
        // assert that with incomplete arguments overriding that the inbreeding coefficient is not calculated
        Assert.assertEquals(((InbreedingCoeff) incompleteArgumentsOverridingGoodArguments.get(0)).annotate(null, inbreedingCoefficientVC, null), Collections.emptyMap());
    }

    //TODO This test is intended to show that the disabling and reenabling an annotation will restore global defaults but unfortunately
    //     this is not the case due to an open issue having to do with instantiating only one annotation object.
    //     See https://github.com/broadinstitute/gatk/issues/3848 for the same issue in read filters.
    @Test (enabled = false)
    public void testDisableDefaultsAndReplaceOnCommandLine() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), Arrays.asList(new InbreedingCoeff(Sets.newSet("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10")),
                        new DepthPerSampleHC()), Collections.singletonList(StandardAnnotation.class))),
                Collections.emptySet());
        List<String> args = Stream.of(StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, StandardAnnotation.class.getSimpleName(), "--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS).collect(Collectors.toList());

        clp.parseArguments(nullMessageStream, args.toArray(new String[args.size()]));
        VariantAnnotatorEngine vae = new VariantAnnotatorEngine(instantiateAnnotations(clp), null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(DepthPerSampleHC.class.getSimpleName())));
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(BaseQualityRankSumTest.class.getSimpleName())));

        // Asserting that the default set InbreedingCoeff has been overridden
        Assert.assertEquals(Double.valueOf((String) ( vae.annotateContext(inbreedingCoefficientVC, new FeatureContext(),null, null, s -> true))
                        .getAttribute(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)),
                -0, 0.001, "InbreedingCoefficientScores");
    }

    @Test
    public void testNoAnnotations(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());
        String[] args = {};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testGetClassForPluginHelpFromReflection(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());
        String[] args = {};
        clp.parseArguments(nullMessageStream, args);

        final GATKAnnotationPluginDescriptor pluginDescriptor = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);

        Assert.assertTrue(pluginDescriptor.getClassForPluginHelp(RMSMappingQuality.class.getSimpleName()).isInstance(new RMSMappingQuality()));
        Assert.assertTrue(pluginDescriptor.getClassForPluginHelp("Foo")==null);
    }

    @Test
    public void testIncludeDefaultExcludeIndividual(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), Collections.singletonList(new Coverage()), null)),
                Collections.emptySet());
        String[] args = {"-AX", Coverage.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testIncludeDefaultGroupExcludeIndividual(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, Collections.singletonList(StandardAnnotation.class))),
                Collections.emptySet());
        String[] args = {"--annotations-to-exclude", Coverage.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        //check that Coverage is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
    }

    @Test
    public void testIncludeUserGroupExcludeIndividual(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, null)),
                Collections.emptySet());
        String[] args = {"--annotations-to-exclude", Coverage.class.getSimpleName(), "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, StandardAnnotation.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
        //check that Coverage is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
    }

    @Test
    public void testAnnotationGroupOverriding(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), null, Collections.singletonList(StandardHCAnnotation.class))),
                Collections.emptySet());
        String[] args = {"-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, StandardAnnotation.class.getSimpleName()};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, null, Collections.emptyList(), false);

        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        //check that Coverage is in and ClippingRankSum is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
        Assert.assertTrue(vae.getInfoAnnotations().stream().anyMatch(a -> a.getClass().getSimpleName().equals(ClippingRankSumTest.class.getSimpleName())));
    }

    @DataProvider
    public Object[][] groupHierarchyArguments(){
        return new Object[][]{
                // we must add --DISABLE_TOOL_DEFAULT_ANNOTATIONS because we pass our annotation objects in to have their groups dynamically discovered
                {new String[]{"--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS}, false, false},
                {new String[]{"--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME,"ParentAnnotationGroup"}, true, true},
                {new String[]{"--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME,"ChildAnnotationGroup"}, true, false},
                {new String[]{"--"+ StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, "-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME,"ParentAnnotationGroup","-"+StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME,"ChildAnnotationGroup"}, true, true},
                {new String[]{}, true, true},};
    }

    // Though current annotation groups aren't implemented in a hierarchical manner, this test enforces that they would behave
    // in an expected manner as well as testing that annotation groups can be dynamically discovered
    @Test (dataProvider = "groupHierarchyArguments")
    public void testGroupHierarchyBehavior(String[] args, boolean includeChild, boolean includeParent){
        testChildAnnotation childAnnotation = new testChildAnnotation();
        testParentAnnotation parentAnnotation = new testParentAnnotation();
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(), Arrays.asList(parentAnnotation, childAnnotation), null)),
                Collections.emptySet());
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);
        Assert.assertEquals(annots.contains(childAnnotation), includeChild);
        Assert.assertEquals(annots.contains(parentAnnotation), includeParent);
    }

    @Test
    public void testOverridingInstancesWithGetInstance() throws InstantiationException, IllegalAccessException {
        GATKAnnotationPluginDescriptor pluginDescriptor = new GATKAnnotationPluginDescriptor(getEmptyArgs(), Collections.singletonList(new testParentAnnotation(true)), null);
        pluginDescriptor.createInstanceForPlugin(testParentAnnotation.class);

        Collection<Annotation> finalAnnotations = pluginDescriptor.getResolvedInstances();
        Assert.assertEquals(finalAnnotations.size(), 1);

        VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Arrays.asList(finalAnnotations.toArray(new Annotation[0])), null, Collections.emptyList(), false);
        VariantContext vc = inbreedingCoefficientVC;
        vc = vae.annotateContext(vc, new FeatureContext(), null, null, a->true);

        Assert.assertEquals(vc.getAttribute("Parent"), null);
    }

    private interface ParentAnnotationGroup extends Annotation { }
    private interface ChildAnnotationGroup extends ParentAnnotationGroup { }
    static class testChildAnnotation extends InfoFieldAnnotation implements ChildAnnotationGroup  {
        public int argument = 5;

        @Override
        public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, ReadLikelihoods<Allele> likelihoods) {
            return Collections.singletonMap("Child",Integer.toString(argument));
        }
        @Override
        public List<String> getKeyNames() {
            return Collections.singletonList("Test");
        }
    }
    static class testParentAnnotation extends InfoFieldAnnotation implements ParentAnnotationGroup  {
        boolean dontAnnotate = false;

        testParentAnnotation() { }

        testParentAnnotation(boolean toAnnotate) {
            this.dontAnnotate = toAnnotate;
        }

        @Override
        public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, ReadLikelihoods<Allele> likelihoods) {
            if (!dontAnnotate) {
                return Collections.singletonMap("Parent","foo");
            } else {
                return Collections.emptyMap();
            }
        }
        @Override
        public List<String> getKeyNames() {
            return Collections.singletonList("Test");
        }
    }

    @Test
    public void testAll(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(),null, null)),
                Collections.emptySet());
        String[] args = {"--use-all-annotations"};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);

        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, dbSNPBinding, features, false);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions(false);
        Assert.assertFalse(vcfAnnotationDescriptions.isEmpty());

        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY)));          //yes DP
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY)));   //yes AC
    }

    @Test
    public void testAllMinusCoverage(){
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(getEmptyArgs(),null, null)),
                Collections.emptySet());
        String[] args = {"--use-all-annotations", "-AX", "Coverage"};
        clp.parseArguments(nullMessageStream, args);
        List<Annotation> annots = instantiateAnnotations(clp);

        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annots, dbSNPBinding, features, false);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions(false);
        Assert.assertFalse(vcfAnnotationDescriptions.isEmpty());

        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY))); //no DP
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY)));  //yes AC
    }
}