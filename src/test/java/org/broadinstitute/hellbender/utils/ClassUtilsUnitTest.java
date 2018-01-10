package org.broadinstitute.hellbender.utils;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotation;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class ClassUtilsUnitTest extends GATKBaseTest {

    //These classes are for testing the methods in ClassUtils.

    private static class Priv{
    }

    public static class PubNestedUtils{
        private PubNestedUtils(){}
    }

    public static class PubNested{
        public PubNested(){}
    }

    public static class PubNestedImplicitConstructor{
    }

    static class NestedWithPublicConstructor{
        public NestedWithPublicConstructor(){}
    }

    static class NestedImplicitConstructor{
    }
    @Test
    public void testCanMakeInstances() {
        Assert.assertFalse(ClassUtils.canMakeInstances(Collections.class));//a util class
        Assert.assertFalse(ClassUtils.canMakeInstances(PubNestedUtils.class));
        Assert.assertFalse(ClassUtils.canMakeInstances(NestedImplicitConstructor.class));

        Assert.assertTrue(ClassUtils.canMakeInstances(PubNested.class));
        Assert.assertTrue(ClassUtils.canMakeInstances(PubNestedImplicitConstructor.class));
        Assert.assertTrue(ClassUtils.canMakeInstances(NestedWithPublicConstructor.class));

        Assert.assertFalse(ClassUtils.canMakeInstances(null));
        Assert.assertFalse(ClassUtils.canMakeInstances(int.class));
        Assert.assertFalse(ClassUtils.canMakeInstances(List.class));

        class Local{}
        Assert.assertFalse(ClassUtils.canMakeInstances(Local.class));

        Assert.assertFalse(ClassUtils.canMakeInstances(Priv.class));

        Assert.assertFalse(ClassUtils.canMakeInstances(AbstractList.class));
        Assert.assertTrue(ClassUtils.canMakeInstances(ArrayList.class));
    }

    @Test
    public void testMakeInstancesOfAnnotations() throws Exception {
        //check that this does not blow up - ie each VariantAnnotation has a non-arg constructor
        ClassUtils.makeInstancesOfSubclasses(VariantAnnotation.class, VariantAnnotation.class.getPackage());
    }

    public abstract static class C{}
    public static class C1 extends C{}
    public static final class C2 extends C{}
    public static final class C11 extends C1{}

    @Test
    public void testMakeInstances() throws Exception {
        final List<? extends C> cs = ClassUtils.makeInstancesOfSubclasses(C.class, C.class.getPackage());
        final Set<String> actualSet = cs.stream().map(o -> o.getClass().getSimpleName()).collect(Collectors.toSet());
        final Set<String> expectedSet = Sets.newHashSet(C1.class.getSimpleName(), C11.class.getSimpleName(), C2.class.getSimpleName());
        Assert.assertEquals(actualSet, expectedSet);
    }

    @Test
    public void testMakeInstancesOfCLP() throws Exception {
        //check that this does not blow up - ie CLPs have a non-arg constructor
        ClassUtils.makeInstancesOfSubclasses(CommandLineProgram.class, ApplyBQSR.class.getPackage());
    }


    interface A{}
    interface A1 extends A{}
    interface B1 extends A1{}

    @Test
    public void testKnownSubinterfaces() throws Exception {
        Assert.assertEquals(new LinkedHashSet<>(ClassUtils.knownSubInterfaceSimpleNames(A.class)), new LinkedHashSet<>(Arrays.asList("A1", "B1")));
    }

    @DataProvider
    private Object[][] provideForTestGetClassesOfType() {

        final List<Class<?>> mapSubclasses  = new ArrayList<>( ClassUtils.knownSubInterfaces(Map.class) );
        final List<Class<?>> listSubclasses = new ArrayList<>( ClassUtils.knownSubInterfaces(List.class) );

        return new Object[][] {
                {
                        Map.class,
                        Stream.concat(mapSubclasses.stream(), listSubclasses.stream()).collect(Collectors.toList()),
                        mapSubclasses
                },
                {
                        List.class,
                        Stream.concat(mapSubclasses.stream(), listSubclasses.stream()).collect(Collectors.toList()),
                        listSubclasses
                },
                {
                        B1.class,
                        new ArrayList<>( ClassUtils.knownSubInterfaces(B1.class) ),
                        new ArrayList<>()
                },
        };
    }

    @Test(dataProvider = "provideForTestGetClassesOfType")
    public void testGetClassesOfType_list( final Class<?> clazz, final List<Class<?>> classesToCheck, final List<Class<?>> expected) {
        Assert.assertEquals( ClassUtils.getClassesOfType(clazz, classesToCheck), expected );
    }
}
