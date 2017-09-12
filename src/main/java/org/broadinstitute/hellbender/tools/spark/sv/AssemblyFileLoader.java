package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.reflect.io.Path;

import java.util.List;
import java.util.WeakHashMap;

/**
 * Created by valentin on 9/11/17.
 */
public class AssemblyFileLoader {

    private final String fastqDir;
    private final String fastqFileFormat;


    private final WeakHashMap<String, List<Template>> contents;

    public AssemblyFileLoader(final String fastqDir, final String fastqFileFormat) {
        this.fastqDir = Utils.nonNull(fastqDir);
        this.fastqFileFormat = Utils.nonNull(fastqFileFormat);
        this.contents = new WeakHashMap<>(10);
    }

    public List<Template> get(final int assemblyNumber) {
        final String assemblyFile = String.format(fastqFileFormat, assemblyNumber);
        final InputStream in = BucketUtils.openFile()

    }
}
