package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.LogEvent;
import org.apache.logging.log4j.core.Logger;
import org.apache.logging.log4j.core.appender.AbstractAppender;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ReferenceConfidenceVariantContextMergerMLEACUnitTest {
    @Test
    public void testStaleMLEAnnotationsAreSkippedBeforeNumericParsing() throws IOException {
        final Path input = Files.createTempFile("mleac-merge-", ".g.vcf");
        try {
            Files.writeString(input, "##fileformat=VCFv4.2\n"
                    + "##contig=<ID=1,length=1000>\n"
                    + "##ALT=<ID=NON_REF,Description=\"Any alternative allele\">\n"
                    + "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele counts\">\n"
                    + "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"MLE allele frequencies\">\n"
                    + "##INFO=<ID=X,Number=1,Type=Integer,Description=\"Other numeric annotation\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods\">\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    + "1\t100\t.\tA\tC,<NON_REF>\t.\tPASS\tMLEAC=1,0;MLEAF=0.500,0.00;X=7\tGT:PL\t0/1:100,0,100,200,200,300\n");
            try (final VCFFileReader reader = new VCFFileReader(input, false)) {
                final VariantContext vc = reader.iterator().next();
                Assert.assertTrue(vc.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY) instanceof List,
                        "The fixture must exercise the decoded Number=A list");
                final VariantAnnotatorEngine engine = new VariantAnnotatorEngine(
                        Collections.emptyList(), null, Collections.emptyList(), false, false);
                final ReferenceConfidenceVariantContextMerger merger =
                        new ReferenceConfidenceVariantContextMerger(engine, reader.getFileHeader());
                final Logger logger = (Logger) LogManager.getLogger(ReferenceConfidenceVariantContextMerger.class);
                final WarningCapture warnings = new WarningCapture();
                warnings.start();
                logger.addAppender(warnings);
                final VariantContext merged;
                try {
                    merged = merger.merge(Collections.singletonList(vc), new SimpleInterval("1", 100, 100),
                            (byte) 'A', true, false);
                } finally {
                    logger.removeAppender(warnings);
                    warnings.stop();
                }
                Assert.assertNotNull(merged);
                Assert.assertFalse(merged.hasAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
                Assert.assertFalse(merged.hasAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
                Assert.assertEquals(merged.getAttributeAsInt("X", -1), 7);
                Assert.assertTrue(warnings.messages.stream().noneMatch(s -> s.contains("Detected invalid annotations")),
                        "Valid stale MLE fields must not produce an invalid-annotation warning: " + warnings.messages);
            }
        } finally {
            Files.deleteIfExists(input);
        }
    }

    private static final class WarningCapture extends AbstractAppender {
        private final List<String> messages = new ArrayList<>();

        private WarningCapture() {
            super("mleac-warning-capture", null, null, false);
        }

        @Override
        public void append(final LogEvent event) {
            if (event.getLevel().isMoreSpecificThan(Level.WARN)) {
                messages.add(event.getMessage().getFormattedMessage());
            }
        }
    }
}
