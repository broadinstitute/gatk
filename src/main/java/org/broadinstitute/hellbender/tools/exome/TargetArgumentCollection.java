package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Set of arguments for tools that work with targets.
 *
 * <p>
 *   The main argument, {@link #targetsFile}, is the targets containing file name itself.
 * </p>
 * <p>
 *   In addition the user can specify several target annotation files: {@link #targetAnnotations}.
 * </p>
 * <p>
 *   If the user does not specify such a file name, the default target file supplier is used to
 *   determine a suitable replacement.
 * </p>
 * <p>
 *   If even this mechanism does not provide a target file name, and the target file is deemed non-optional,
 *   the a {@link UserException.BadArgumentValue} is thrown.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetArgumentCollection {

    public static final String TARGET_FILE_SHORT_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME;
    public static final String TARGET_FILE_LONG_NAME = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME;
    public static final String TARGET_ANNOTATION_FILES_SHORT_NAME = "A";
    public static final String TARGET_ANNOTATION_FILES_FULL_NAME = "targetAnnotations";

    @Argument(
            doc = "Targets listing file",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_LONG_NAME,
            optional = true
    )
    protected File targetsFile;

    @Argument(
            doc = "Target annotations",
            shortName = TARGET_ANNOTATION_FILES_SHORT_NAME,
            fullName = TARGET_ANNOTATION_FILES_FULL_NAME,
            optional = true
    )
    protected List<File> targetAnnotations = new ArrayList<>();

    /**
     * Holds a reference to a function that provides a default target file in case the user
     * does not specify one explicitly using argument {@link #targetsFile}.
     */
    private final Supplier<File> defaultTargetsFileSupplier;

    /**
     * Creates a new target argument collection.
     * <p>
     *    The one argument for this constructor is a reference to the default target file supplier.
     * </p>
     *
     * @param defaultTargetsFileSupplier the fail-over targets file name supplier in case no target
     *                                    file name was indicated explicitly.
     * @throws IllegalArgumentException if {@code defaultTargetsFileSupplier} is {@code null}.
     */
    public TargetArgumentCollection(final Supplier<File> defaultTargetsFileSupplier) {
        this.defaultTargetsFileSupplier = Utils.nonNull(defaultTargetsFileSupplier, "the default target file supplier cannot be null");
    }

    /**
     * Creates a new target argument collection without a default target file supplier.
     */
    public TargetArgumentCollection() {
        this(() -> null);
    }

    /**
     * Reads in a target collection from the target file provided by the user.
     * <p>
     *     The only argument for this function, {@code optional}, determines whether a target collection
     *     is at all required by the enclosing tool and what happens when this is missing:
     * </p>
     * <ul>
     * <li><p>
     *     If {@code optional == false} a {@link UserException} is thrown prompting the user to use the
     *     {@link #targetsFile} argument in order to specify one explicitly.
     * </p></li>
     * <li><p>
     *     If {@code optional == true} a {@code null} is returned.
     * </p></li>
     * </lu>
     *
     * @param optional whether the target collection is optional or not.
     * @return may be {@code null}, if {@code optional} is {@code true} and no target file could be
     *      determined.
     * @throws UserException.BadArgumentValue if {@code optional} is {@code false} and no target file
     *   was provided or cannot be resolved otherwise.
     */
    public TargetCollection<Target> readTargetCollection(final boolean optional) {
        final File resolveFile = getTargetsFile();
        if (resolveFile != null) {
            if (targetAnnotations.isEmpty()) {
                return readTargetCollection(resolveFile);
            } else {
                return combineTargetAnnotations(resolveFile, targetAnnotations);
            }
        } else {
            if (optional) {
                return null;
            } else {
                throw new UserException.BadArgumentValue(TARGET_FILE_SHORT_NAME, "No target file was specified");
            }
        }
    }

    /**
     * Reads a target collection from a file and adds all the annotations in additional
     * annotation files.
     * @param targetFile the target collection file.
     * @param targetAnnotations additional annotation files.
     * @return never {@code null}.
     */
    private TargetCollection<Target> combineTargetAnnotations(final File targetFile,
                                                              final List<File> targetAnnotations) {
        final TargetCollection<Target> targets = readTargetCollection(targetFile);

        final Map<Target, Map<TargetAnnotation, String>> annotationsByTarget =
                composeTargetToAnnotationMap(targetAnnotations, targets);

        final List<Target> fullyAnnotatedTargets = targets.targets().stream()
                .map(originalTarget -> new Target(originalTarget.getName(),
                                originalTarget.getInterval(),
                                new TargetAnnotationCollection(annotationsByTarget.get(originalTarget))))
                .collect(Collectors.toList());

        return new HashedListTargetCollection<Target>(fullyAnnotatedTargets) {

            @Override
            public String name(final Target target) {
                return target.getName();
            }

            @Override
            public SimpleInterval location(final Target target) {
                return target.getInterval();
            }
        };
    }

    /**
     * Composes a Map from target to its annotations.
     * <p>
     *     The output map will contain exactly one entry per target in
     *     the input target-collection {@code targets}. The value map
     *     will contain all the annotations (keys) for each target along with their values.
     *     Such maps may be empty if there is no annotation values for that target.
     * </p>
     *
     * @param targetAnnotations files containing target annotations.
     * @param targets the targets to annotate.
     * @return never {@code null.}
     */
    private Map<Target, Map<TargetAnnotation, String>> composeTargetToAnnotationMap(final List<File> targetAnnotations,
                                                                                    final TargetCollection<Target> targets) {

        final Map<Target, Map<TargetAnnotation, String>> annotationsByTarget =
                targets.targets().stream().collect(Collectors.toMap(target -> target, TargetArgumentCollection::targetToTargetAnnotationMap));

        for (final File annotationFile : targetAnnotations) {
            try (final TargetTableReader reader = new TargetTableReader(annotationFile)) {
                reader.stream()
                        .filter(annotationsByTarget::containsKey)
                        .forEach(annotatedTarget -> {
                            for (final TargetAnnotation annotation : annotatedTarget.getAnnotations().annotationSet()) {
                                final String value = annotatedTarget.getAnnotations().get(annotation);
                                annotationsByTarget.get(annotatedTarget).put(annotation, value);
                            }
                        });
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(annotationFile, ex);
            }
        }
        return annotationsByTarget;
    }

    /**
     * Composes a pair from the input target to a mutable map with is annotation values.
     * @param target the input target.
     * @return never {@code null}. Value maps are never {@code null} either, but they might be empty
     * if the input target does not have any annotations.
     */
    private static Map<TargetAnnotation, String> targetToTargetAnnotationMap(final Target target) {

        final Map<TargetAnnotation, String> result = new EnumMap<>(TargetAnnotation.class);
        for (final TargetAnnotation annotation : target.getAnnotations().annotationSet()) {
            result.put(annotation, target.getAnnotations().get(annotation));
        }
        return result;
    }

    /**
     * Resolves the name of the target file given the explicit argument {@link #targetsFile} and the
     * default target file supplier {@link #defaultTargetsFileSupplier}.
     *
     * @return {@code null} when neither approach resolves into a target file name. The file returned is not
     * guaranteed to exist or to be a readable regular file.
     */
    public File getTargetsFile() {
        return targetsFile != null ? targetsFile
                : (defaultTargetsFileSupplier != null ) ? defaultTargetsFileSupplier.get() : null;
    }

    /**
     * Reads the content of a file into a targets collection.
     * @param file the file to read.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem when reading the file
     *         provided.
     */
    public static TargetCollection<Target> readTargetCollection(final File file) {
        try (final TargetTableReader reader = new TargetTableReader(file)) {
            return new HashedListTargetCollection<Target>(Utils.nonNull(reader.stream().collect(Collectors.toList()), "the input feature list cannot be null")) {
                @Override
                public String name(final Target target) {
                    return Utils.nonNull(target,"the input target cannot be null").getName();
                }

                @Override
                public SimpleInterval location(final Target target) {
                    return Utils.nonNull(target, "the input target cannot be null").getInterval();
                }
            };
        } catch (final IOException | UncheckedIOException ex) {
            throw new UserException.CouldNotReadInputFile(file, ex.getMessage());
        }
    }
}
