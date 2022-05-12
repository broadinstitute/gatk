package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONObject;

import java.io.*;
import java.nio.file.Path;

public abstract class RampBase {

    protected final Logger logger = LogManager.getLogger(getClass());

    // type of ramp enum - determines data flow direction
    public enum Type {
        OffRamp,
        OnRamp
    }

    // local vars
    final protected Type        type;
    final protected File        file;
    protected JSONObject        info;

    public RampBase(final String filename, final Type type) throws IOException  {
        this.type = type;
        this.file = new File(filename);
        logger.info("opening ramp. file: " + file + ", type: " + type);
    }

    protected String getLocFilenameSuffix(final Locatable loc) {
        return String.format("%s-%d-%d", loc.getContig(), loc.getStart(), loc.getEnd());
    }

    public void close() throws IOException {
    }

    public Type getType() {
        return type;
    }

    protected Path getBamIndexPath(Path bamPath) {
        return (new File(bamPath.toFile().getAbsolutePath().replace(".bam", ".bai"))).toPath();
    }

    protected void copyStreamToPath(final InputStream is, final Path outPath) throws IOException {
        final OutputStream        os = new FileOutputStream(outPath.toFile());
        IOUtils.copy(is, os);
        os.close();
    }

    protected String getReadSuppName(final GATKRead read) {
        return getReadSuppName(read.getName(), read.isSupplementaryAlignment());
    }

    protected String getReadSuppName(final String name, final boolean supp) {
        return name + (supp ? ",1" : ",0");
    }


}
