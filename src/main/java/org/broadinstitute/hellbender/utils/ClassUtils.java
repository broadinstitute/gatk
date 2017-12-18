package org.broadinstitute.hellbender.utils;

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Utilities for dealing with reflection.
 */
public final class ClassUtils {
    private ClassUtils(){}

    /**
     * Returns true iff we can make instances of this class.
     * Note that this will return false if the class does not have any public constructors.
     */
    public static boolean canMakeInstances(final Class<?> clazz) {
        return clazz != null &&
                !clazz.isPrimitive()  &&
                !clazz.isSynthetic()  &&
                !clazz.isInterface()  &&
                !clazz.isLocalClass() &&
                !Modifier.isPrivate(clazz.getModifiers()) &&
                !Modifier.isAbstract(clazz.getModifiers()) &&
                clazz.getConstructors().length != 0;
    }

    /**
     * Finds and creates objects of all concrete subclasses of the given class in the package.
     * The public no-arg constructor is called to create the objects.
     *
     * GATKException is thrown if creation of any object fails.
     * @param clazz class to be instantiated
     * @param pack package in which the class will be searched for
     */
    @SuppressWarnings("unchecked")
    public static <T> List<T> makeInstancesOfSubclasses(final Class<? extends T> clazz, final Package pack){
        Utils.nonNull(clazz, "class");
        Utils.nonNull(pack, "package");
        final ClassFinder finder = new ClassFinder();
        finder.find(pack.getName(), clazz);
        final Set<Class<?>> classes = finder.getClasses();

        final List<T> results = new ArrayList<>(classes.size());

        for (final Class<?> found: classes){
            if (canMakeInstances(found)){
                try {
                    results.add((T)found.newInstance());
                } catch (InstantiationException | IllegalAccessException e) {
                    throw new GATKException("Problem making an instance of " + found + " Do check that the class has a non-arg constructor", e);
                }
            }
        }
        return results;
    }

    /**
     * Finds sub-interfaces of the given interface (in the same package) and returns their simple names.
     */
    public static List<String> knownSubInterfaceSimpleNames(final Class<?> iface) {
        Utils.nonNull(iface);
        Utils.validateArg(iface.isInterface(), iface + " is not an interface");
        return knownSubInterfaces(iface).stream().map(c ->c.getSimpleName()).collect(Collectors.toList());
    }

    /**
     * Finds all subinterfaces of the given interface (in the same package).
     */
    public static Set<Class<?>> knownSubInterfaces(final Class<?> iface) {
        final ClassFinder finder = new ClassFinder();
        finder.find(iface.getPackage().getName(), iface);
        return finder.getClasses().stream().filter(cl -> !cl.equals(iface) && cl.isInterface()).collect(Collectors.toSet());
    }

    /**
     * Gets a list of classes that are either the same as, or a subclass/subinterface of a parent target class.
     * @param targetClass Parent {@link Class} for which to check for inheritance.
     * @param classesToSearch Classes to check for inheritance against {@code targetClass}.
     * @return {@link List} of classes from {@code classesToSearch} that inherit from {@code targetClass}.
     */
    public static List<Class<?>> getClassesOfType(final Class<?> targetClass, final List<Class<?>> classesToSearch) {

        final List<Class<?>> classList = new ArrayList<>();

        for ( final Class<?> clazz : classesToSearch ) {
            if ( targetClass.isAssignableFrom(clazz) ) {
                classList.add( clazz );
            }
        }

        return classList;
    }
}
