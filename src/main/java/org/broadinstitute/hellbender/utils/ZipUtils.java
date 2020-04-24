package org.broadinstitute.hellbender.utils;

import com.google.cloud.storage.Acl;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URI;
import java.util.Arrays;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

public class ZipUtils {

    public static final void unzip(final String path, final File dest, final String ... files) {
        final Predicate<File> selected = composeSelectPredicate(files);
        final File file = stageLocallyIfNeeded(path);
        try (final ZipFile zipFile = new ZipFile(file)) {
            Utils.stream(zipFile.entries())
                    .filter(entry -> selected.test(new File(entry.getName())))
                    .forEach(entry -> {
                        final File entryDest = new File(dest, entry.getName());
                        if (entry.isDirectory()) {
                            entryDest.mkdirs();
                        } else {
                            final File parent = entryDest.getParentFile();
                            if (!parent.exists()) {
                                parent.mkdirs();
                            }
                            try {
                                IOUtils.copy(zipFile.getInputStream(entry), new FileOutputStream(entryDest));
                            } catch (final IOException ex) {
                                throw new GATKException("could not dicompress temporary file " + entryDest);
                            }
                        }
                    });
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(file, ex);
        }
    }

    private static File stageLocallyIfNeeded(String path) {
        final File file;
        if (BucketUtils.isFileUrl(path)) {
            file = new File(URI.create(path).getPath());
            if (!file.exists()) {
                throw new UserException.CouldNotReadInputFile(file);
            } else if (!file.isFile()) {
                throw new UserException.CouldNotReadInputFile("" + file + " is not a regular file (e.g. directory)");
            } else if (!file.canRead()) {
                throw new UserException.CouldNotReadInputFile("" + file + " is not readable");
            }
        } else {
            try {
                file = File.createTempFile("download", ".zip");
                file.deleteOnExit();
                BucketUtils.copyFile(path, file.toString());
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(path, ex);
            }
        }
        return file;
    }

    private static Predicate<File> composeSelectPredicate(final String[] files) {
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
