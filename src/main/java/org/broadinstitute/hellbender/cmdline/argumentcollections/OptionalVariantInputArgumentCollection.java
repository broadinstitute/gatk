package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * An argument collection for use with tools that accept zero or more input files containing VariantContext records
 * (eg., VCF files).
 */
public final class OptionalVariantInputArgumentCollection extends VariantInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              doc = "Variants file", optional = true)
    public List<FeatureInput<VariantContext>> variantFiles;
    public List<String> variantFileNames;

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    @Override
    public List<Path> getVariantPaths() {
      ArrayList<Path> ret = new ArrayList<>();
      for (String fn : variantFileNames) {
        ret.add(IOUtils.getPath(fn));
      }
      return ret;
    }

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line
     */
    @Override
    public List<File> getVariantFiles() {
      ArrayList<File> ret = new ArrayList<>();
      for (String fn : variantFileNames) {
        ret.add(new File(fn));
      }
      return ret;
    }
}
