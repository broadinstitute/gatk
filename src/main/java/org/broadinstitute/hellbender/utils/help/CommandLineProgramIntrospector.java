package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLinePluginProvider;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.Utils;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Utility to query definitions of program arguments given their names.
 */
public final class CommandLineProgramIntrospector {

    private final static Map<Class<?>, CommandLineProgramIntrospector> introspectors = new HashMap<>();
    private static final BiMap<String, Class<?>> toolToClass = HashBiMap.create();

    /*
     * TODO this code was copied from Picard code base,
     * TODO we need to expose this functionality in Picard and then pull it from its the dependency
     * TODO instead.
     */
    static {

        final List<String> packageList = Arrays.asList("org.broadinstitute.hellbender", "picard");

        final ClassFinder classFinder = new ClassFinder();
        for (final String pkg : packageList) {
            classFinder.find(pkg, picard.cmdline.CommandLineProgram.class);
            classFinder.find(pkg, CommandLineProgram.class);
        }
        final Set<Class<?>> toCheck = classFinder.getClasses();
        for (final Class<?> clazz : toCheck) {
        		if (!ClassUtils.canMakeInstances(clazz)) {
        			continue;
        		}
        		final CommandLineProgramProperties property = clazz.getAnnotation(CommandLineProgramProperties.class);
            if (property == null) {
            		continue;
            }
            toolToClass.put(clazz.getSimpleName(), clazz);
        }
    }

    private final CommandLineArgumentParser parser;
    private final Map<String, CommandLineArgumentParser.ArgumentDefinition> definitionsByName;

    /**
     * Create an introspector instance given the tool class.
     * @param clazz the target tool class.
     * @throws IllegalArgumentException if we cannot create an instance of the tool 
     * class using the parameter-less constructor.
     */
    private CommandLineProgramIntrospector(final Class<?> clazz) {
        this.definitionsByName = new HashMap<>();
        try {
            final Object program = clazz.newInstance();
            final CommandLineArgumentParser parser;
            if (program instanceof CommandLinePluginProvider) {
                final List<? extends CommandLinePluginDescriptor<?>> pluginDescriptors = ((CommandLinePluginProvider) program).getPluginDescriptors();
                parser = new CommandLineArgumentParser(program, pluginDescriptors, Collections.emptySet());
            } else {
                parser = new CommandLineArgumentParser(program);
            }
            this.parser = parser;
        } catch (final InstantiationException | IllegalAccessException e) {
            throw new IllegalArgumentException("could not create parser for " + clazz.getName(), e);
        }
        for (final CommandLineArgumentParser.ArgumentDefinition definition : parser.getArgumentDefinitions()) {
            for (final String name : definition.getNames()) {
                definitionsByName.put(name, definition);
            }
        }
    }

    /**
     * Returns the introspector for a tool given its name.
     * @param name the target tool name.
     * @return null if there is no tool with such a name.
     */
    public static CommandLineProgramIntrospector of(final String name) {
        final Class<?> clazz = toolToClass.get(name);
        return clazz == null 
        		? null
        		: introspectors.computeIfAbsent(clazz, CommandLineProgramIntrospector::new);
    }

    /**
     * Returns the introspector for a tool given its class.
     * @param clazz the target tool class.
     * @throws IllegalArgumentException if {@code clazz} is {@code null}.
     * @return never {@code null}.
     */
    public static CommandLineProgramIntrospector of(final Class<?> clazz) {
    		return clazz == null 
    				? null
    				: introspectors.computeIfAbsent(clazz, CommandLineProgramIntrospector::new);
    }

    /**
     * Returns the {@link CommandLineArgumentParser.ArgumentDefinition ArgumentDefinition} given either
     *  the short or full user argument name for the target tool.
     * @param name either the short or full user argument name of interest.
     * @return {@code null} iff there is no argument with such a name.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     */
    public CommandLineArgumentParser.ArgumentDefinition getArgumentDefinition(final String name) {
        Utils.nonNull(name);
        return definitionsByName.get(name);
    }
}
