/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.diffengine;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;


/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 1:09 PM
 *
 * Class implementing diffnode reader for VCF
 */
public class VCFDiffableReader implements DiffableReader {
    private static Logger logger = Logger.getLogger(VCFDiffableReader.class);

    @Override
    public String getName() { return "VCF"; }

    @Override
    public DiffElement readFromFile(File file, int maxElementsToRead) {
        DiffNode root = DiffNode.rooted(file.getName());
        try {
            // read the version line from the file
            BufferedReader br = new BufferedReader(new FileReader(file));
            final String version = br.readLine();
            root.add("VERSION", version);
            br.close();

            final VCFCodec vcfCodec = new VCFCodec();
            vcfCodec.disableOnTheFlyModifications(); // must be read as state is stored in reader itself

            FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(file.getAbsolutePath(), vcfCodec, false);
            VCFHeader header = (VCFHeader)reader.getHeader();
            for ( VCFHeaderLine headerLine : header.getMetaDataInInputOrder() ) {
                String key = headerLine.getKey();
                if ( headerLine instanceof VCFIDHeaderLine)
                    key += "_" + ((VCFIDHeaderLine) headerLine).getID();
                if ( root.hasElement(key) )
                    logger.warn("Skipping duplicate header line: file=" + file + " line=" + headerLine.toString());
                else
                    root.add(key, headerLine.toString());
            }

            int count = 0, nRecordsAtPos = 1;
            String prevName = "";
            Iterator<VariantContext> it = reader.iterator();
            while ( it.hasNext() ) {
                VariantContext vc = it.next();
                String name = vc.getChr() + ":" + vc.getStart();
                if ( name.equals(prevName) ) {
                    name += "_" + ++nRecordsAtPos;
                } else {
                    prevName = name;
                }
                DiffNode vcRoot = DiffNode.empty(name, root);

                // add fields
                vcRoot.add("CHROM", vc.getChr());
                vcRoot.add("POS", vc.getStart());
                vcRoot.add("ID", vc.getID());
                vcRoot.add("REF", vc.getReference());
                vcRoot.add("ALT", vc.getAlternateAlleles());
                vcRoot.add("QUAL", vc.hasLog10PError() ? vc.getLog10PError() * -10 : VCFConstants.MISSING_VALUE_v4);
                vcRoot.add("FILTER", ! vc.filtersWereApplied() // needs null to differentiate between PASS and .
                        ? VCFConstants.MISSING_VALUE_v4
                        : ( vc.getFilters().isEmpty() ? VCFConstants.PASSES_FILTERS_v4 : vc.getFilters()) );

                // add info fields
                for (Map.Entry<String, Object> attribute : vc.getAttributes().entrySet()) {
                    if ( ! attribute.getKey().startsWith("_") )
                        vcRoot.add(attribute.getKey(), attribute.getValue());
                }

                for (Genotype g : vc.getGenotypes() ) {
                    DiffNode gRoot = DiffNode.empty(g.getSampleName(), vcRoot);
                    gRoot.add("GT", g.getGenotypeString());
                    if ( g.hasGQ() ) gRoot.add("GQ", g.getGQ() );
                    if ( g.hasDP() ) gRoot.add("DP", g.getDP() );
                    if ( g.hasAD() ) gRoot.add("AD", Utils.join(",", g.getAD()));
                    if ( g.hasPL() ) gRoot.add("PL", Utils.join(",", g.getPL()));
                    if ( g.getFilters() != null ) gRoot.add("FT", g.getFilters());

                    for (Map.Entry<String, Object> attribute : g.getExtendedAttributes().entrySet()) {
                        if ( ! attribute.getKey().startsWith("_") )
                            gRoot.add(attribute.getKey(), attribute.getValue());
                    }

                    vcRoot.add(gRoot);
                }

                root.add(vcRoot);
                count += vcRoot.size();
                if ( count > maxElementsToRead && maxElementsToRead != -1)
                    break;
            }

            reader.close();
        } catch ( IOException e ) {
            return null;
        }

        return root.getBinding();
    }

    @Override
    public boolean canRead(File file) {
        return AbstractVCFCodec.canDecodeFile(file.getPath(), VCFCodec.VCF4_MAGIC_HEADER);
    }
}
