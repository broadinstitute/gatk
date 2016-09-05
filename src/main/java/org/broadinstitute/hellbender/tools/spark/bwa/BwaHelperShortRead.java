package org.broadinstitute.hellbender.tools.spark.bwa;

import com.github.lindenb.jbwa.jni.ShortRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * The purpose of this data structure is to enable keeping track of which ShortReads correspond to which GATKReads.
 * The JNI layer ignores the read field.
 */
public final class BwaHelperShortRead extends ShortRead{
    private final GATKRead read; //may be null

    public BwaHelperShortRead(String name, byte[] bases, byte[] baseQualities, GATKRead read) {
        super(name, bases, baseQualities);
        this.read = read;
    }

    public GATKRead getRead(){ return read; }
}
