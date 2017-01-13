package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * An argument collection for use with tools that accept one or more input files containing VariantContext records
 * (eg., VCF files), and require at least one such input.
 */
public final class RequiredVariantInputArgumentCollection extends VariantInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more files containing variants", optional = false)
    public List<FeatureInput<VariantContext>> variantFiles;

    public List<String> variantFilesNames;

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    @Override
    public List<Path> getReadPaths() {
        ArrayList<Path> ret = new ArrayList<>();
        for (String fn : variantFilesNames) {
            ret.add(IOUtils.getPath(fn));
        }
        return ret;
    }

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line
     */
    @Override
    public List<File> getReadFiles() {
        ArrayList<File> ret = new ArrayList<>();
        for (String fn : variantFilesNames) {
            ret.add(new File(fn));
        }
        return ret;
    }

    /**
     * Get the list of BAM/SAM/CRAM filenames specified at the command line
     */
    public List<String> getVariantFilesNames() {
        return new ArrayList<>(variantFilesNames);
    }
}