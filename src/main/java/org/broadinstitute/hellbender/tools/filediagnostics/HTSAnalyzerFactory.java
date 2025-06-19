package org.broadinstitute.hellbender.tools.filediagnostics;

import htsjdk.samtools.util.FileExtensions;
import org.broadinstitute.hellbender.engine.GATKPath;

/**
 * Class for creating an analyzer based on an alignment file type.
 */
public class HTSAnalyzerFactory {

    public static HTSAnalyzer getFileAnalyzer(final GATKPath inputPath, final GATKPath outputPath, final long countLimit) {
        System.out.println(inputPath.getRawInputString());
        if (inputPath.isCram()) {
            return new CRAMAnalyzer(inputPath, outputPath, countLimit);
        } else if (inputPath.hasExtension(FileExtensions.CRAM_INDEX)) {
            return new CRAIAnalyzer(inputPath, outputPath);
        } else if (inputPath.hasExtension(FileExtensions.BAI_INDEX)) {
            return new BAIAnalyzer(inputPath, outputPath);
        } else {
            throw new RuntimeException("Unsupported diagnostic file type: " + inputPath.getRawInputString());
        }
    }
}
