package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Checks on Caller argument collection cloning.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class StandardCallerArgumentCollectionUnitTest {

    public final static List<Class<? extends StandardCallerArgumentCollection>> COLLECTION_CLASSES;

    static {
        COLLECTION_CLASSES = new ArrayList<>();
        COLLECTION_CLASSES.add(StandardCallerArgumentCollection.class);
        COLLECTION_CLASSES.add(UnifiedArgumentCollection.class);
    }

    @Test(dataProvider="collectionClasses")
    public void testParameterLessInitialization(final Class<? extends StandardCallerArgumentCollection> clazz) {
           try {
               clazz.newInstance();
           } catch (final Exception ex) {
               Assert.fail();
           }
    }

    @DataProvider(name="collectionClasses")
    public Object[][] collectionClasses() {
        final Object[][] result = new Object[COLLECTION_CLASSES.size()][];
        for (int i = 0; i < COLLECTION_CLASSES.size(); i++)
            result[i] = new Object[] { COLLECTION_CLASSES.get(i) };
        return result;
    }

    @DataProvider(name="collectionClassPairs")
    public Object[][] collectionClassPairs() {
        final Object[][] collectionClasses = collectionClasses();
        final Object[][] result = new Object[collectionClasses.length * collectionClasses.length][];
        for (int i = 0; i < collectionClasses.length; i++)
            for (int j = 0; j < collectionClasses.length; j++)
                result[i * collectionClasses.length + j] = new Object[] { collectionClasses[i][0], collectionClasses[j][0]};
        return result;
    }

    public <T extends StandardCallerArgumentCollection> T randomArgumentCollection(final Class<T> clazz) throws IllegalAccessException, InstantiationException {
        final T result = clazz.newInstance();
        final Random rnd = Utils.getRandomGenerator();
        for (final Field field : clazz.getFields()) {
            final int fieldModifiers = field.getModifiers();
            if (!Modifier.isPublic(fieldModifiers))
                continue;
            final Class<?> fieldType = mapWrappersToPrimitives(field.getType());
            final Object value;
            if (fieldType.isPrimitive()) {
                if (fieldType == Integer.TYPE)
                    value = rnd.nextInt(100) - 50;
                else if (fieldType == Double.TYPE)
                    value = rnd.nextDouble() * 10;
                else if (fieldType == Float.TYPE)
                    value = rnd.nextFloat() * 10;
                else if (fieldType == Boolean.TYPE)
                    value = rnd.nextBoolean();
                else if (fieldType == Byte.TYPE)
                    value = (byte) rnd.nextInt(256);
                else if (fieldType == Character.TYPE)
                    value = (char) rnd.nextInt(256);
                else if (fieldType == Long.TYPE)
                    value = rnd.nextLong();
                else if (fieldType == Short.TYPE)
                    value = (short) rnd.nextLong();
                else
                    throw new IllegalStateException("unknown primitive type!!!!");
            } else if (fieldType == String.class) {
                value = "" + rnd.nextLong();
            } else if (fieldType.isEnum()) {
                value = fieldType.getEnumConstants()[rnd.nextInt(fieldType.getEnumConstants().length)];
            } else { // Cannot handle other types in a general way.
                value = null;
            }
            field.set(result,value);
        }
        return result;
    }

    private Class<?> mapWrappersToPrimitives(Class<?> type) {
        if (type == Boolean.class)
            return Boolean.TYPE;
        else if (type == Integer.class)
            return Integer.TYPE;
        else if (type == Double.class)
            return Double.TYPE;
        else if (type == Float.class)
            return Float.TYPE;
        else if (type == Character.class)
            return Character.TYPE;
        else if (type == Short.class)
            return Short.TYPE;
        else if (type == Byte.class)
            return Byte.TYPE;
        else if (type == Long.class)
            return Long.TYPE;
        else
            return type;
    }
}
