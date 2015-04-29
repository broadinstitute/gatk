package org.broadinstitute.hellbender.utils.diffengine;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.Iterator;
import java.util.Map;


/**
 * Class implementing diffnode reader for VCF
 */
final class VCFDiffableReader implements DiffableReader {
    private static Logger logger = LogManager.getLogger(VCFDiffableReader.class);

    @Override
    public String getName() { return "VCF"; }

    @Override
    public DiffElement readFromFile(final File file, final int maxElementsToRead) throws IOException {
        final DiffNode root = DiffNode.rooted(file.getName());
        // read the version line from the file
        final BufferedReader br = new BufferedReader(new FileReader(file));
        final String version = br.readLine();
        root.add("VERSION", version);
        br.close();

        final VCFCodec vcfCodec = new VCFCodec();
        vcfCodec.disableOnTheFlyModifications(); // must be read as state is stored in reader itself

        final FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), vcfCodec, false);
        final VCFHeader header = (VCFHeader)reader.getHeader();
        for ( final VCFHeaderLine headerLine : header.getMetaDataInInputOrder() ) {
            String key = headerLine.getKey();
            if ( headerLine instanceof VCFIDHeaderLine) {
                key += "_" + ((VCFIDHeaderLine) headerLine).getID();
            }
            if ( root.hasElement(key) ) {
                logger.warn("Skipping duplicate header line: file=" + file + " line=" + headerLine.toString());
            } else {
                root.add(key, headerLine.toString());
            }
        }

        int count = 0, nRecordsAtPos = 1;
        String prevName = "";
        final Iterator<VariantContext> it = reader.iterator();
        while ( it.hasNext() ) {
            final VariantContext vc = it.next();
            String name = vc.getContig() + ":" + vc.getStart();
            if ( name.equals(prevName) ) {
                name += "_" + ++nRecordsAtPos;
            } else {
                prevName = name;
            }
            final DiffNode vcRoot = DiffNode.empty(name, root);

            // add fields
            vcRoot.add("CHROM", vc.getContig());
            vcRoot.add("POS", vc.getStart());
            vcRoot.add("ID", vc.getID());
            vcRoot.add("REF", vc.getReference());
            vcRoot.add("ALT", vc.getAlternateAlleles());
            vcRoot.add("QUAL", vc.hasLog10PError() ? vc.getLog10PError() * -10 : VCFConstants.MISSING_VALUE_v4);
            vcRoot.add("FILTER", ! vc.filtersWereApplied() // needs null to differentiate between PASS and .
                    ? VCFConstants.MISSING_VALUE_v4
                    : ( vc.getFilters().isEmpty() ? VCFConstants.PASSES_FILTERS_v4 : vc.getFilters()) );

            // add info fields
            for (final Map.Entry<String, Object> attribute : vc.getAttributes().entrySet()) {
                if ( ! attribute.getKey().startsWith("_") ) {
                    vcRoot.add(attribute.getKey(), attribute.getValue());
                }
            }

            for (final Genotype g : vc.getGenotypes() ) {
                final DiffNode gRoot = DiffNode.empty(g.getSampleName(), vcRoot);
                gRoot.add("GT", g.getGenotypeString());
                if ( g.hasGQ() ) {
                    gRoot.add("GQ", g.getGQ());
                }
                if ( g.hasDP() ) {
                    gRoot.add("DP", g.getDP());
                }
                if ( g.hasAD() ) {
                    gRoot.add("AD", Utils.join(",", g.getAD()));
                }
                if ( g.hasPL() ) {
                    gRoot.add("PL", Utils.join(",", g.getPL()));
                }
                if ( g.getFilters() != null ) {
                    gRoot.add("FT", g.getFilters());
                }

                for (final Map.Entry<String, Object> attribute : g.getExtendedAttributes().entrySet()) {
                    if ( ! attribute.getKey().startsWith("_") ) {
                        gRoot.add(attribute.getKey(), attribute.getValue());
                    }
                }

                vcRoot.add(gRoot);
            }

            root.add(vcRoot);
            count += vcRoot.size();
            if ( count > maxElementsToRead && maxElementsToRead != -1) {
                break;
            }
        }

        reader.close();
        return root.getBinding();
    }

    @Override
    public boolean canRead(final File file) {
        return AbstractVCFCodec.canDecodeFile(file.getPath(), VCFCodec.VCF4_MAGIC_HEADER);
    }
}
