package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Set of arguments for tools that work with targets.
 *
 * <p>
 *   The main argument is the targets containing file name itself.
 * </p>
 * <p>
 *   If the user does not specify such a file name, the fail-over target file supplier is used to
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

    public static final String TARGET_FILE_SHORT_NAME = "T";
    public static final String TARGET_FILE_FULL_NAME = "targets";

    @Argument(
            doc = "Targets listing file",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_FULL_NAME,
            optional = true
    )
    protected File targetsFile;

    /**
     * Holds a reference to a function that provides a default target file in case the user
     * does not specify one explicitly using argument {@link #targetsFile}.
     */
    private final Supplier<File> defaultTargetFileSupplier;

    /**
     * Creates a new target argument collection.
     * <p>
     *    The one argument for this constructor is a reference to the default target file supplier.
     * </p>
     *
     * @param defaultTargetFileSupplier the fail-over targets file name supplier in case no target
     *                                    file name was indicated explicitly.
     * @throws IllegalArgumentException if {@code defaultTargetFileSupplier} is {@code null}.
     */
    protected TargetArgumentCollection(final Supplier<File> defaultTargetFileSupplier) {
        this.defaultTargetFileSupplier = Utils.nonNull(defaultTargetFileSupplier, "the default target file supplier cannot be null");
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
    protected TargetCollection<Target> readTargetCollection(final boolean optional) {
        final File resolveFile = getTargetsFile();
        if (resolveFile == null) {
            if (optional) {
                return null;
            } else {
                throw new UserException.BadArgumentValue(TARGET_FILE_SHORT_NAME,
                        "No target file was specified");
            }
        } else {
            return readTargetCollection(resolveFile);
        }
    }


    /**
     * Resolves the name of the target file given the explicit argument {@link #targetsFile} and the
     * default target file supplier {@link #defaultTargetFileSupplier}.
     *
     * @return {@code null} when neither approach resolves into a target file name. The file returned is not
     * guaranteed to exist or been a readable regular file.
     */
    public File getTargetsFile() {
        return targetsFile != null ? targetsFile
                : (defaultTargetFileSupplier != null ) ? defaultTargetFileSupplier.get()
                : null;
    }

    /**
     * Reads the content of a file into a targets collection.
     * @param file the file to read.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem when reading the file
     *         provided.
     */
    private static TargetCollection<Target> readTargetCollection(final File file) {
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
