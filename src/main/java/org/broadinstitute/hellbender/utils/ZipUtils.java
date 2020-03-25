package org.broadinstitute.hellbender.utils;

import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.*;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.zip.*;

/**
 * Utility class to zip and unzip files.
 */
public class ZipUtils {

    /**
     * Unzips a file in a provided path into a destination directory.
     * @param path the location of the zip file to unzip.
     * @param dest the destination directory.
     * @param files files to unzip, by omission is the whole zip file content.
     */
    public static void unzip(final GATKPath path, final File dest, final String ... files) {
        Utils.nonNull(path, "the path provided cannot be null");
        Utils.nonNull(dest, "the destination directory provided cannot be null");
        Utils.validateArg((dest.isDirectory() || dest.mkdirs()) && dest.canWrite(), "could not create destination directory");
        final Predicate<File> selected = composeSelectPredicate(files);
        try (final ZipInputStream in = new ZipInputStream(path.getInputStream())) {
            ZipEntry entry;
            while ((entry = in.getNextEntry()) != null) {
                if (selected.test(new File(entry.getName()))) {
                    final File entryDest = new File(dest, entry.getName());
                    if (entry.isDirectory()) {
                        entryDest.mkdirs();
                    } else {
                        final File parent = entryDest.getParentFile();
                        if (!parent.exists()) {
                            parent.mkdirs();
                        }
                        try (final FileOutputStream fos = new FileOutputStream(entryDest) ) {
                            IOUtils.copy(in, fos);
                        } catch (final IOException e) {
                            throw new GATKException("problems unzipping entry", e);
                        }
                    }
                }
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(path.toString(), ex);
        }
    }

    /**
     * Creates a zip file give a source file/directory to zip and the final destination zip file path.\
     * <p>
     *     Directories are zipped recursively.
     * </p>
     * @param source the source folder or file.
     * @param dest the destination path.
     * @param files if not empty, this array indicate what files needs to be zipped. If zero length, all
     *               files will be zipped.
     */
    public static void zip(final File source, final GATKPath dest, final String ... files) {
        Utils.nonNull(source, "the destination directory provided cannot be null");
        Utils.nonNull(dest, "the destination cannot be null");
        if (source.isFile()) {
            if (files != null && files.length != 0) {
                throw new IllegalArgumentException("if the source is a regular file, you cannot specify a files to zip array");
            }
            uncheckedZip(source.getParentFile(), dest, composeSelectPredicate(source.toString()));
        } else {
            final Predicate<File> selection = composeSelectPredicate(files);
            uncheckedZip(source, dest, selection);
        }
    }

    /**
     * Performs the unzip action assuming inputs are correct.
     */
    private static void uncheckedZip(final File sourceRoot, final GATKPath dest, final Predicate<File> selection) {
        try (final ZipOutputStream out = new ZipOutputStream(dest.getOutputStream())) {
            final Deque<File> pending = new ArrayDeque<>();
            final String prefix = sourceRoot.toString() + File.separator;
            pending.addAll(Arrays.asList(sourceRoot.listFiles()));
            while (!pending.isEmpty()) {
                final File next = pending.removeFirst();
                if (next.isFile()) {
                    if (selection.test(next)) {
                        final String relative  = next.toString().replace(prefix, "");
                        out.putNextEntry(new ZipEntry(relative));
                        try (final FileInputStream nextIs = new FileInputStream(next)) {
                            IOUtils.copy(nextIs, out);
                        }
                        out.closeEntry();
                    }
                } else if (next.isDirectory()) {
                    final File[] children = next.listFiles();
                    for (int i = children.length - 1; i >= 0; i--) {
                        pending.addFirst(children[i]);
                    }
                }
            }
        } catch (final IOException ex) {
            throw new GATKException("errors trying to create zip file from " + sourceRoot +
                     " to " + dest, ex);
        }
    }

    /**
     * Returns a predicate that returns {@code true} for entries present in the input array of {@code files}
     * @param files the target files.
     * @return never {@code null}.
     */
    private static Predicate<File> composeSelectPredicate(final String ... files) {
        if (files == null || files.length == 0) {
            return _x -> true;
        } else {
            final Set<File> selected = Arrays.stream(files)
                    .map(name -> new File(name).toPath().normalize().toFile())
                    .collect(Collectors.toSet());
            return f -> selected.contains(f) || selected.contains(f.toPath().normalize().toFile());
        }
    }
}
