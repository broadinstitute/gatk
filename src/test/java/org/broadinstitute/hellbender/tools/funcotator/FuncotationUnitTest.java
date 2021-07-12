package org.broadinstitute.hellbender.tools.funcotator;

import com.esotericsoftware.kryo.Kryo;
import com.google.common.collect.Sets;
import htsjdk.utils.ClassFinder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.GATKRegistrator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.test.FuncotatorTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedHashSet;

public class FuncotationUnitTest extends GATKBaseTest {

    private static final String FUNCOTATION_ROOT_PACKAGE_FOR_SEARCH = "org.broadinstitute.hellbender.tools.funcotator";

    /**
     * If the serialization fails, address the error message as seen in
     *  {@link org.broadinstitute.hellbender.engine.spark.GATKRegistrator#registerFuncotationMapDependencies(Kryo)}
     * We need to test that all concrete subclasses of Funcotation can be serialized in Kryo.
     * If this is not the case (or classes are missing from the test).  Then Fail.
     * Note that this relies on all concrete implementations of Funcotation implement equals().
     *
     */
    @Test
    public void testAllFuncotationConcreteClassCanSerialize() {

        // Add future implementations of Funcotation to the map of class to instance value.
        final LinkedHashMap<Class<?>, Funcotation> funcotationToInstanceMap = FuncotatorUtils.createLinkedHashMapFromLists(

                // Note: If you need to include another serialization test, add to the two lists below (the lists must have corresponding entries).
                //  The first list is the concrete class instance (E.g. FooFuncotation.class) and the second list is an example instance.

                Arrays.asList(
                        GencodeFuncotation.class,
                        TableFuncotation.class
                ),
                Arrays.asList(
                        // The values do not matter that much here.
                        FuncotatorTestUtils.createDummyGencodeFuncotation("TEST_TX", new VariantContextBuilder().chr("1").start(10000).stop(10000).alleles("C", "A").make()),
                        FuncotatorTestUtils.createDummyTableFuncotation()
                )
        );

        /* Test that we are testing all concrete classes that extend Funcotation. */
        final ClassFinder classFinder = new ClassFinder();
        // Query class finder for classes
        classFinder.find(FUNCOTATION_ROOT_PACKAGE_FOR_SEARCH, Funcotation.class);

        final Set<Class<?>> concreteFuncotationClasses = classFinder.getConcreteClasses();
        Assert.assertTrue(concreteFuncotationClasses.size() > 0, "No concrete funcotation classes found, which indicates that this test is misconfigured.");

        final Set<Class<?>> implementedButNotFound = Sets.difference(funcotationToInstanceMap.keySet(), concreteFuncotationClasses).immutableCopy();
        final Set<Class<?>> notImplementedButFound = Sets.difference(concreteFuncotationClasses, funcotationToInstanceMap.keySet()).immutableCopy();

        Assert.assertTrue(implementedButNotFound.size() == 0, "The following Funcotation concrete classes were in the test, but not actually found in " + FUNCOTATION_ROOT_PACKAGE_FOR_SEARCH + ": " +
          implementedButNotFound.stream().map(Class::getName).collect(Collectors.joining(", ")));

        Assert.assertTrue(notImplementedButFound.size() == 0, "The following Funcotation concrete classes were not included in this test, but found in " + FUNCOTATION_ROOT_PACKAGE_FOR_SEARCH + ": " +
          implementedButNotFound.stream().map(Class::getName).collect(Collectors.joining(", ")) + ".  Please add to the method FuncotationUnitTest::testAllFuncotationConcreteClassCanSerialize");

        Assert.assertEquals(funcotationToInstanceMap.keySet(), concreteFuncotationClasses);

        /* Test that the given instances can be serialized as a round trip.  */
        final Consumer<Kryo> registerFuncotationMapForKryo = GATKRegistrator::registerFuncotationMapDependencies;
        for (final Class<?> concreteFuncotationClass : funcotationToInstanceMap.keySet()) {
            final Funcotation testInstance = funcotationToInstanceMap.get(concreteFuncotationClass);
            final Funcotation testInstanceCopy = FuncotatorTestUtils.assertRoundTripInKryo(testInstance, concreteFuncotationClass,
                    Collections.singletonList(registerFuncotationMapForKryo));
            Assert.assertFalse(testInstance == testInstanceCopy);
            Assert.assertEquals(testInstance, testInstanceCopy);
        }
    }


}

