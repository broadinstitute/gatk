//package org.broadinstitute.hellbender.utils.help;
//
//import htsjdk.samtools.util.Iso8601Date;
//import org.apache.commons.lang3.tuple.ImmutablePair;
//import org.apache.commons.lang3.tuple.Pair;
//
//import org.broadinstitute.barclay.argparser.*;
//import org.broadinstitute.barclay.help.DocWorkUnit;
//import org.broadinstitute.barclay.help.HelpDoclet;
//import org.broadinstitute.barclay.help.TemplateProperties;
//import org.broadinstitute.barclay.help.WDLWorkUnitHandler;
//import org.broadinstitute.hellbender.engine.FeatureInput;
//import org.broadinstitute.hellbender.engine.GATKPath;
//import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
//import org.broadinstitute.hellbender.exceptions.GATKException;
//import picard.illumina.parser.ReadStructure;
//
//import java.util.*;
//
//// Note: WDL Gen doesn't handle arguments that accept tagged argument values
//
///**
// * The GATK WDL work unit handler. Its main task is to convert the types for all arguments for a given work
// * unit (tool) from Java types to WDL-compatible types by updating the freemarker map with the transformed types.
// *
// * NOTE: Methods in this class are intended to be called by Gradle/Javadoc only, and should not be called
// * by methods that are used by the GATK runtime, as this class assumes a dependency on com.sun.javadoc classes
// * which may not be present.
// */
//public class GATKWDLWorkUnitHandler extends WDLWorkUnitHandler {
//
//    private final static String GATK_FREEMARKER_TEMPLATE_NAME = "wdlToolTemplate.wdl.ftl";
//
//    // This must be kept in sync with the value used in build.gradle, where the file is created
//    private final static String dummyWDLTestFileName = "dummyWDLTestFile";
//
//    private final static String WDL_TEST_VALUE_PROPERTY = "testValue";
//
//    // Map of Java argument types that the WDL generator knows how to convert to a WDL type, along with the
//    // corresponding string substitution that needs to be run on the (Barclay-generated) string that describes
//    // the type. From a purely string perspective, some of these transforms are no-ops in that no actual
//    // conversion is required because the type names are identical in Java and WDL (i.e, File->File or
//    // String->String), but they're included here for completeness, and to document the allowed type transitions.
//    private final static Map<Class<?>, ImmutablePair<String, String>> javaToWDLTypeMap =
//            new HashMap<Class<?>, ImmutablePair<String, String>>() {
//                private static final long serialVersionUID = 1L;
//                {
//                    // GATK-specific File Types
//                    put(GATKPath.class, new ImmutablePair<>(GATKPath.class.getSimpleName(), "File"));
//                    // FeatureInputs require special handling to account for the generic type param(s)
//                    put(FeatureInput.class, new ImmutablePair<>(FeatureInput.class.getSimpleName(), "File"));
//
//                    put(Iso8601Date.class, new ImmutablePair<>(Iso8601Date.class.getSimpleName(), "String"));
//                    put(Date.class, new ImmutablePair<>(Date.class.getSimpleName(), "String"));
//                    put(ReadStructure.class, new ImmutablePair<>(ReadStructure.class.getSimpleName(), "String"));
//                }
//        };
//
//    public GATKWDLWorkUnitHandler(final HelpDoclet doclet) {
//        super(doclet);
//    }
//
//    /**
//     * @param workUnit the DocWorkUnit object being processed
//     * @return the name of a the freemarker template to be used for the class being documented.
//     * Must reside in the folder passed to the Barclay Doclet via the "-settings-dir" parameter to
//     * Javadoc.
//     */
//    @Override
//    public String getTemplateName(final DocWorkUnit workUnit) { return GATK_FREEMARKER_TEMPLATE_NAME; }
//
//    /**
//     * Return the flat filename (no paths) that the handler would like the Doclet to
//     * write out the documentation for workUnit.
//     * @param workUnit
//     * @return the name of the destination file to which documentation output will be written
//     */
//    @Override
//    public String getDestinationFilename(final DocWorkUnit workUnit) {
//        return workUnit.getClazz().getSimpleName() + ".wdl";
//    }
//
//    /**
//     * Returns the JSON output file name.
//     */
//    @Override
//    public String getJSONFilename(final DocWorkUnit workUnit) {
//        return workUnit.getClazz().getSimpleName() + "Inputs.json";
//    }
//
//    /**
//     * Add the named argument {@code argDef}to the property map if applicable.
//     * @param currentWorkUnit current work unit
//     * @param args the freemarker arg map
//     * @param argDef the arg to add
//     */
//    @Override
//    @SuppressWarnings("unchecked")
//    protected void processNamedArgument(
//            final DocWorkUnit currentWorkUnit,
//            final Map<String, List<Map<String, Object>>> args,
//            final NamedArgumentDefinition argDef)
//    {
//        // for WDL gen, we don't want the special args such as --help or --version to show up in the
//        // WDL or JSON input files
//        if (!argDef.getUnderlyingField().getDeclaringClass().equals(SpecialArgumentsCollection.class)) {
//            super.processNamedArgument(currentWorkUnit, args, argDef);
//
//            // now extract any newlines out of the summary since the summary appears in a WDL line comment
//            final List<Map<String, Object>> argMapList = args.get("all");
//            argMapList.stream().forEach(
//                    m -> {
//                        final String actualArgName = (String) m.get(TemplateProperties.WDL_ARGUMENT_ACTUAL_NAME);
//                        if (actualArgName != null && actualArgName.equals("--" + argDef.getLongName())) {
//                            // the summary string must be suitable to embed in a quoted string in the param_meta section,
//                            // so remove any quotes and newlines
//                            final String newSummary = ((String) m.get(TemplateProperties.ARGUMENT_SUMMARY)).replaceAll("[\"\'\n]", "");
//                            m.put(TemplateProperties.ARGUMENT_SUMMARY, newSummary);
//                        }
//                    });
//        }
//    }
//
//    @Override
//    protected String processNamedArgument(
//            final Map<String, Object> argBindings,
//            final NamedArgumentDefinition argDef,
//            final String fieldCommentText) {
//        final String argCategory = super.processNamedArgument(argBindings, argDef, fieldCommentText);
//        argBindings.put(WDL_TEST_VALUE_PROPERTY,
//                testValueAsJSON(
//                        argDef.getLongName(),
//                        argDef,
//                        (String) argBindings.get(TemplateProperties.ARGUMENT_TYPE),
//                        (String) argBindings.get(TemplateProperties.ARGUMENT_DEFAULT_VALUE))
//        );
//
//        return argCategory;
//    }
//
//    @Override
//    protected void processPositionalArguments(
//            final CommandLineArgumentParser clp,
//            final Map<String, List<Map<String, Object>>> args) {
//        super.processPositionalArguments(clp, args);
//        final List<Map<String, Object>> positionalArgsList = args.get(TemplateProperties.ARGUMENTS_POSITIONAL);
//        if (positionalArgsList != null && !positionalArgsList.isEmpty()) {
//            final Map<String, Object> positionalArgs = args.get(TemplateProperties.ARGUMENTS_POSITIONAL).get(0);
//            positionalArgs.put(WDL_TEST_VALUE_PROPERTY,
//                    testValueAsJSON(
//                            WDLWorkUnitHandler.POSITIONAL_ARGS,
//                            clp.getPositionalArgumentDefinition(),
//                            (String) positionalArgs.get(TemplateProperties.ARGUMENT_TYPE),
//                            (String) positionalArgs.get(TemplateProperties.ARGUMENT_DEFAULT_VALUE))
//            );
//        }
//    }
//
//    /**
//     * Return a test input value for use in the WDL validation test inputs.
//     *
//     * If an option has WDL type File, then we need to provide the name of an actual file that exists so cromwell
//     * can localize it: "dummyWDLTestFileName" is a file that will be created by the test task.
//     *
//     * @param longName the long name for this arg
//     * @param argDef the ArgumentDefinition for this arg
//     * @param wdlType the wdl type for which an input value is needed
//     * @param defaultWDLValue the default value for the argument for which an input value is required
//     * @return a test input value that is either the default value, or the name of an actual test file
//     * that will exist at test execution time
//     */
//    protected String testValueAsJSON(
//            final String longName,
//            final ArgumentDefinition argDef,
//            final String wdlType,
//            final String defaultWDLValue) {
//        final Argument argumentAnnotation = argDef.getUnderlyingField().getAnnotation(Argument.class);
//        final PositionalArguments positionalAnnotation = argDef.getUnderlyingField().getAnnotation(PositionalArguments.class);
//        final boolean isRequired =
//                (argumentAnnotation != null && !argumentAnnotation.optional())
//                || positionalAnnotation != null;
//        final String dummyWDLTestFile = ((GATKWDLDoclet) getDoclet()).getBuildDir() + "/" + dummyWDLTestFileName;
//
//        // Hack to resolve the mutex argument in GATKSparkTool; otherwise all spark tools will fail since all
//        // mutex args will have a value
//        if (longName.equals(GATKSparkTool.OUTPUT_SHARD_DIR_LONG_NAME)) {
//            return "null";
//        }
//
//        // first check the wdl type; for File we always want to use the name of the dummy file that is created by
//        // the test process that is used to ensure localization is handled correctly
//        if (wdlType.equals("File")) {
//            if (isRequired) {
//                return "\"" + dummyWDLTestFile + "\"";
//            } else {
//                return "null";
//            }
//        } else if (wdlType.equals("Array[File]")) {
//            return String.format("[\"%s\", \"%s\"]", dummyWDLTestFile, dummyWDLTestFile);
//        }
//
//        // for other (non-File) types, use the default value and arg def to synthesize a value
//        if (defaultWDLValue.equals("null") || defaultWDLValue.equals("\"\"") || defaultWDLValue.equals("[]")) {
//            if (isRequired) {
//                if (wdlType.startsWith("Array")) {
//                    // we use two values to accommodate tools that take positional args, such as CompareSams that
//                    // require more than one arg
//                    return String.format("[\"%s\", \"%s\"]", getDefaultValueForType(argDef), getDefaultValueForType(argDef));
//                } else {
//                    return "\"" + getDefaultValueForType(argDef) + "\"";
//                }
//            } else {
//                return "null";
//            }
//        } else if (defaultWDLValue.startsWith("[")) {
//            // for required arrays, we need to provide SOME value in the array
//            if (isRequired) {
//                return defaultWDLValue;
//            } else {
//                return "null";
//            }
//        } else if (wdlType.equals("Float")) {
//            if (defaultWDLValue.equals("Infinity") || defaultWDLValue.equals("Nan")) {
//                // JSON does not recognize "Infinity" or "Nan" as valid float values (!), so we
//                // need to treat them as String values
//                return "\"" + defaultWDLValue + "\"";
//            }
//        }
//        return defaultWDLValue;
//    }
//
//    /**
//     * Return a value that will satisfy the constructor for argDef's underlying class.
//     * @param argDef
//     * @return a value that will satisfy the constructor for argDef's underlying class.
//     */
//    protected String getDefaultValueForType(final ArgumentDefinition argDef) {
//        final Class<?> clazz = argDef.getUnderlyingFieldClass();
//        if (clazz.isEnum()) {
//            // any enum constant for this enum will do
//            return argDef.getUnderlyingFieldClass().getEnumConstants()[0].toString();
//        } else {
//            // get a value that is formatted appropriately by calling toString() on the value
//            // on the underlying object if its not null, otherwise return a synthetic string
//            // suitable for the underlying type
//            final Object fieldDefaultValue = argDef.getArgumentValue();
//            if (fieldDefaultValue != null && !fieldDefaultValue.toString().equals("[]")) {
//                return fieldDefaultValue.toString();
//            } else if (clazz.isAssignableFrom(Boolean.class) ||
//                    clazz.isAssignableFrom(boolean.class)) {
//                return "true";
//            } else if (clazz.isAssignableFrom(Float.class) ||
//                    clazz.isAssignableFrom(float.class) ||
//                    clazz.isAssignableFrom(Double.class) ||
//                    clazz.isAssignableFrom(double.class)) {
//                return "0.0";
//            } else if (clazz.isAssignableFrom(Integer.class) ||
//                    clazz.isAssignableFrom(int.class) ||
//                    clazz.isAssignableFrom(Long.class) ||
//                    clazz.isAssignableFrom(long.class)) {
//                return "1";
//            } else {
//                // give up and return a String with a synthetic value
//                return "syntheticTestValue";
//            }
//        }
//    }
//
//    /**
//     * Given a Java class representing the underlying field  type of an argument, and a human readable doc type,
//     * convert the docType to a WDL type.
//     *
//     * @param workflowOutput the WorkflowOutput associated with the instance of argumentClass, if any
//     * @param argumentClass the Class for the underlying field of the argument being converted
//     * @param docType a string representing the human readable type assigned by the Barclay doc system
//     * @param sourceContext a String describing the context for this argument, used for error reporting
//     * @return the docType string transformed to the corresponding WDL type
//     */
//    @Override
//    protected String convertJavaTypeToWDLType(
//            final WorkflowOutput workflowOutput,
//            final Class<?> argumentClass,
//            final String docType,
//            final String sourceContext) {
//        String convertedWDLType;
//        if (FeatureInput.class.isAssignableFrom(argumentClass)) {
//            if (!docType.contains(FeatureInput.class.getSimpleName())) {
//                throw new GATKException(
//                        String.format(
//                                "Don't know how to convert Java type %s in %s to a corresponding WDL type. " +
//                                        "The WDL generator type converter code must be updated to support this Java type.",
//                                argumentClass,
//                                sourceContext));
//            }
//            final Pair<String, String> typeConversionPair = transformToWDLType(argumentClass);
//            convertedWDLType = docType.replaceFirst("FeatureInput\\[[a-zA-Z0-9?]+\\]", typeConversionPair.getValue());
//
//            // finally, if this type is for an arg that is a WorkflowOutput that is a workflow output, and its type
//            // is file, we need to use a different type (String) as the input type for this arg to prevent the workflow
//            // manager from attempting to localize the (non-existent) output file when localizing inputs
//            return transformWorkflowOutputTypeToInputType(workflowOutput, convertedWDLType);
//        }
//        return super.convertJavaTypeToWDLType(workflowOutput, argumentClass, docType, sourceContext);
//    }
//
//    /**
//     * Given an argument class, return a String pair representing the string that should be replaced (the Java type),
//     * and the string to substitute (the corresponding WDL type), i.e., for an argument with type Java Integer.class,
//     * return the Pair ("Integer", "Int") to convert from the Java type to the corresponding WDL type.
//     * @param argumentClass Class of the argument being converter
//     * @return a String pair representing the original and replacement type text, or null if no conversion is available
//     */
//    @Override
//    protected Pair<String, String> transformToWDLType(final Class<?> argumentClass) {
//        Pair<String, String> conversion = javaToWDLTypeMap.get(argumentClass);
//        return conversion == null ? super.transformToWDLType(argumentClass) : conversion;
//    }
//
//    /**
//     * Given a Java collection class, return a String pair representing the string that should be replaced (the Java type),
//     * and the string to substitute (the corresponding WDL type), i.e., for an argument with type Java List.class,
//     * return the Pair ("List", "Array") to convert from the Java type to the corresponding WDL collection type.
//     * @param argumentCollectionClass collection Class of the argument being converter
//     * @return a String pair representing the original and replacement type text, or null if no conversion is available
//     */
//    @Override
//    protected Pair<String, String> transformToWDLCollectionType(final Class<?> argumentCollectionClass) {
//        // required for Picard LiftoverVcf
//        return argumentCollectionClass.equals(Collection.class) ?
//                new ImmutablePair<>("Collection", "Array"):
//                super.transformToWDLCollectionType(argumentCollectionClass);
//    }
//
//    /**
//     * Add any custom freemarker bindings discovered via custom javadoc tags. Subclasses can override this to
//     * provide additional custom bindings.
//     *
//     * @param currentWorkUnit the work unit for the feature being documented
//     */
//    @Override
//    protected void addCustomBindings(final DocWorkUnit currentWorkUnit) {
//        super.addCustomBindings(currentWorkUnit);
//
//        // Picard tools use the summary line for the long overview section, so extract that
//        // from Picard tools only, and put it in the freemarker map.
//        Class<?> toolClass = currentWorkUnit.getClazz();
//        if (picard.cmdline.CommandLineProgram.class.isAssignableFrom(toolClass)) {
//            final CommandLineProgramProperties clpProperties = currentWorkUnit.getCommandLineProperties();
//            currentWorkUnit.setProperty("picardsummary", clpProperties.summary());
//        }
//
//        // add the buildDir as a property so it can be accessed by the test inputs JSON file
//        currentWorkUnit.setProperty("buildDir", ((GATKWDLDoclet)getDoclet()).getBuildDir());
//    }
//
//}
