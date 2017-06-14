package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep reads that do NOT contain one or more kmers from a set of SVKmerShorts
 */
public class ContainsKmerReadFilter extends ReadFilter {

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private static final long serialVersionUID = 1L;
    private static PSKmerCollection kmerLib = null;
    private final int kSize, kmerCountThreshold;

    public ContainsKmerReadFilter(final String kmerLibPath, final int kmerCountThreshold) {
        this.kmerCountThreshold = kmerCountThreshold;
        if (kmerLib == null) {
            synchronized (ContainsKmerReadFilter.class) {
                if (kmerLib == null) {
                    kmerLib = PSKmerUtils.readKmerFilter(kmerLibPath);
                    logger.info("Loaded static kmer filter with false positive probability " + kmerLib.getFalsePositiveProbability());
                }
            }
        }
        kSize = kmerLib.kmerSize();
    }

    @Override
    public boolean test(final GATKRead read) {
        final SVKmerizer kmers = new SVKmerizer(read.getBases(), kSize, 1, new SVKmerShort(kSize));
        int numKmersFound = 0;
        while (kmers.hasNext()) {
            if (kmerLib.contains(((SVKmerShort)kmers.next()))) {
                if (++numKmersFound >= kmerCountThreshold) {
                    return false;
                }
            }
        }
        return true;
    }

    //Static variables can't be garbage collected until the object is unloaded
    public static void closeKmerLib() {
        kmerLib = null;
    }
}
