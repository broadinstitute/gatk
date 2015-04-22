package org.broadinstitute.hellbender.cmdline;


import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.net.URLClassLoader;
import java.net.URLDecoder;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * Utility class that can scan for classes in the classpath and find all the ones
 * annotated with a particular annotation.
 *
 * @author Tim Fennell
 */
public final class ClassFinder {
    private final Set<Class<?>> classes = new HashSet<>();
    private final ClassLoader loader;
    private Class<?> parentType;
    // If not null, only look for classes in this jar
    private String jarPath = null;

    private static final Log log = Log.getInstance(ClassFinder.class);

    public ClassFinder() {
        loader = Thread.currentThread().getContextClassLoader();
    }

    public ClassFinder(final ClassLoader loader) {
        this.loader = loader;
    }

    public ClassFinder(final File jarFile) throws IOException {
        // The class loader must have the context in order to load dependent classes when loading a class,
        // but the jarPath is remembered so that the iteration over the classpath skips anything other than
        // the jarPath.
        jarPath = jarFile.getCanonicalPath();
        final URL[] urls = {new URL("file", "", jarPath)};
        loader = new URLClassLoader(urls, Thread.currentThread().getContextClassLoader());
    }

    /** Convert a filename to a class name by removing '.class' and converting '/'s to '.'s. */
    public String toClassName(final String filename) {
        return filename.substring(0, filename.lastIndexOf(".class"))
                .replace('/', '.').replace('\\', '.');
    }

    /**
     * Scans the classpath for classes within the specified package and sub-packages that
     * extend the parentType. This method can be called repeatedly
     * with different packages. Classes are accumulated internally and
     * can be accessed by calling {@link #getClasses()}.
     */
    public void find(String packageName, final Class<?> parentType) {
        this.parentType = parentType;
        packageName = packageName.replace('.', '/');
        final Enumeration<URL> urls;

        try {
            urls = loader.getResources(packageName);
        }
        catch (IOException ioe) {
            log.warn("Could not read package: " + packageName, ioe);
            return;
        }

        while (urls.hasMoreElements()) {
            try {
                String urlPath = urls.nextElement().getFile();
                urlPath = URLDecoder.decode(urlPath, "UTF-8");
                if ( urlPath.startsWith("file:") ) {
                    urlPath = urlPath.substring(5);
                }
                if (urlPath.indexOf('!') > 0) {
                    urlPath = urlPath.substring(0, urlPath.indexOf('!'));
                }
                if (jarPath != null && !jarPath.equals(urlPath)) {
                    continue;
                }

                //Log.info("Looking for classes in location: " + urlPath);
                final File file = new File(urlPath);
                if ( file.isDirectory() ) {
                    scanDir(file, packageName);
                }
                else {
                    scanJar(file, packageName);
                }
            }
            catch (IOException ioe) {
                log.warn("could not read entries", ioe);
            }
        }
    }

    /**
     * Scans the entries in a ZIP/JAR file for classes under the parent package.
     * @param file the jar file to be scanned
     * @param packagePath the top level package to start from
     */
    protected void scanJar(final File file, final String packagePath) throws IOException {
        final ZipFile zip = new ZipFile(file);
        final Enumeration<? extends ZipEntry> entries = zip.entries();
        while ( entries.hasMoreElements() ) {
            final ZipEntry entry = entries.nextElement();
            final String name = entry.getName();
            if (name.startsWith(packagePath)) {
                handleItem(name);
            }
        }
    }

    /**
     * Scans a directory on the filesystem for classes.
     * @param file the directory or file to examine
     * @param path the package path acculmulated so far (e.g. edu/mit/broad)
     */
    protected void scanDir(final File file, final String path) {
        for ( final File child: file.listFiles() ) {
            final String newPath = (path==null ? child.getName() : path + '/' + child.getName() );
            if ( child.isDirectory() ) {
                scanDir(child, newPath);
            }
            else {
                handleItem(newPath);
            }
        }
    }

    /**
     * Checks an item to see if it is a class and is annotated with the specified
     * annotation.  If so, adds it to the set, otherwise ignores it.
     * @param name the path equivelant to the package + class/item name
     */
    protected void handleItem(final String name) {
        if (name.endsWith(".class")) {
            final String classname = toClassName(name);

            try {
                final Class<?> type = loader.loadClass(classname);
                if (parentType.isAssignableFrom(type)) {
                    this.classes.add(type);
                }
            }
            catch (Throwable t) {
                log.debug("could not load class: " + classname, t);
            }
        }
    }

    /** Fetches the set of classes discovered so far. */
    public Set<Class<?>> getClasses() {
        return this.classes;
    }

    /**
     * Fetches the set of classes discovered so far, subsetted down to concrete (non-abstract/interface) classes only
     *
     * @return subset of classes discovered so far including only concrete (non-abstract/interface) classes
     */
    public Set<Class<?>> getConcreteClasses() {
        Set<Class<?>> concreteClassSet = new HashSet<>();

        for ( Class<?> clazz : classes ) {
            if ( isConcrete(clazz) ) {
                concreteClassSet.add(clazz);
            }
        }

        return concreteClassSet;
    }

    /**
     * Determines whether or not the specified class is concrete (ie., non-abstract and non-interface)
     *
     * @param clazz class to check
     * @return true if the class is neither abstract nor an interface, otherwise false
     */
    public static boolean isConcrete( final Class<?> clazz ) {
        return ! Modifier.isAbstract(clazz.getModifiers()) && ! Modifier.isInterface(clazz.getModifiers());
    }
}