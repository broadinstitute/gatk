package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.*;
import java.nio.file.Path;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public abstract class OffRampBase extends RampBase {

    final protected ZipOutputStream outputZip;

    public OffRampBase(final String filename) throws IOException {
        super(filename, Type.OffRamp);

        // open zip for writing
        this.file.getAbsoluteFile().getParentFile().mkdirs();
        outputZip = new ZipOutputStream(new FileOutputStream(this.file));

        // create info object
        info = new JSONObject();
        info.put("regions", new JSONArray());
    }

    public void close() throws IOException {

        // add info
        addEntry(null,"info.json", info.toString(2).getBytes());

        // close file
        outputZip.close();

        super.close();
    }

    protected void addHaplotypes(final Locatable loc, final String name, final List<Haplotype> value) throws IOException {

        // build text representation
        final ByteArrayOutputStream   os = new ByteArrayOutputStream();
        final PrintWriter             pw = new PrintWriter(os);
        pw.println("contig,start,end,ref,cigar,bases,score,alignmentStartHapwrtRef");
        for ( Haplotype haplotype : value ) {
            final Locatable        hapLoc = haplotype.getGenomeLocation();
            pw.println(String.format("%s,%d,%d,%d,%s,%s,%s,%d",
                    hapLoc.getContig(), hapLoc.getStart(), hapLoc.getEnd(),
                    haplotype.isReference() ? 1 : 0,
                    haplotype.getCigar().toString(),
                    new String(haplotype.getBases()),
                    Double.toString(haplotype.getScore()),
                    haplotype.getAlignmentStartHapwrtRef()
                    )
            );
        }
        pw.close();

        // add entry
        addEntry(loc, name, os.toByteArray());
    }

    protected void add(final Locatable loc, final String name, final Object value) throws IOException {

        // this is a useless default
        addEntry(loc, name, value.toString().getBytes());
    }

    protected void addInfo(final JSONObject info, final String name, final Object obj) {

        JSONObject      parent = info;
        final String[]  toks =  name.split("\\.");
        for ( int tokIndex = 0 ; tokIndex < toks.length - 1 ; tokIndex++ ) {
            if (!parent.has(toks[tokIndex])) {
                parent.put(toks[tokIndex], new JSONObject());
            }
            parent = parent.getJSONObject(toks[tokIndex]);
        }
        parent.put(toks[toks.length - 1], obj);
    }

    protected void addEntry(final Locatable loc, final String name, final byte[] bytes) throws IOException {

        final String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
        final ZipEntry    e = new ZipEntry(prefix + name);
        outputZip.putNextEntry(e);
        outputZip.write(bytes);
        outputZip.closeEntry();
    }

    protected void addEntry(final Locatable loc, final String name, final Path path) throws IOException {

        final String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
        final ZipEntry    e = new ZipEntry(prefix + name);
        outputZip.putNextEntry(e);
        final InputStream is = new FileInputStream(path.toFile());
        IOUtils.copy(is, outputZip);
        is.close();
        outputZip.closeEntry();
    }
}
