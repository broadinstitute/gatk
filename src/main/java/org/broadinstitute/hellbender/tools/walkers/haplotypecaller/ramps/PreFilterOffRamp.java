package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONObject;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.stream.Collectors;

public class PreFilterOffRamp extends OffRampBase {

    static final boolean WIDE_READS = false;

    public PreFilterOffRamp(final String filename) throws IOException {
        super(filename);
    }

    public synchronized  void add(final Locatable loc, final String name, final AlleleLikelihoods<GATKRead, Haplotype> value,
                                  final AssemblyRegion regionForGenotyping, final AssemblyRegion region) throws IOException {

        // add region
        final JSONObject regionObj = new JSONObject();
        regionObj.put("contig", loc.getContig());
        regionObj.put("start", loc.getStart());
        regionObj.put("end", loc.getEnd());

        // add global info
        addInfo(regionObj , name + ".haplotypeCount", value.numberOfAlleles());
        addInfo(regionObj, name + ".readCount", value.evidenceCount());
        addInfo(regionObj, name + ".r4g", regionForGenotyping.getSpan());
        addInfo(regionObj, name + ".r4gp", regionForGenotyping.getPaddedSpan());

        // add haplotypes, contigs
        addHaplotypes(loc,name + ".haplotypes", value.alleles());

        // loop on samples
        for ( int sampleIndex = 0 ; sampleIndex < value.numberOfSamples() ; sampleIndex++ ) {

            // establish context
            final String                                 baseName = name + ".samples." + value.getSample(sampleIndex);
            final LikelihoodMatrix<GATKRead, Haplotype> sampleMatrix = value.sampleMatrix(sampleIndex);

            // add to info
            final JSONObject      info2 = new JSONObject();
            info2.put("readCount", sampleMatrix.evidenceCount());
            addInfo(regionObj, baseName, info2);

            // write reads
            final String      msg = regionForGenotyping.toString() + " " + region.toString();
            addReads(loc,baseName + ".reads", sampleMatrix.evidence(), regionForGenotyping.getReads(), region.getReads(), msg);

            // write matrix itself
            addSampleMatrix(loc, baseName + ".matrix", sampleMatrix);
        }

        // add to regions
        info.getJSONArray("regions").put(regionObj);
    }

    private void addSampleMatrix(final Locatable loc, final String name, final LikelihoodMatrix<GATKRead, Haplotype> value) throws IOException {

        // build text representation
        final ByteArrayOutputStream os = new ByteArrayOutputStream();
        final PrintWriter pw = new PrintWriter(os);

        // header line
        pw.print("read,supp");
        for ( int n = 0 ; n < value.numberOfAlleles() ; n++ )
            pw.print(",h" + n);
        pw.println("");

        // walk the lines
        for ( int readIndex  = 0 ; readIndex < value.evidenceCount() ; readIndex++ ) {
            final GATKRead        read = value.getEvidence(readIndex);
            pw.print(read.getName() + "," + (read.isSupplementaryAlignment() ? 1 : 0));
            for ( int hapIndex = 0 ; hapIndex < value.numberOfAlleles() ; hapIndex++ ) {
                //pw.print("," + String.format(MATRIX_FLOAT_FORMAT, value.get(hapIndex, readIndex)));
                pw.print("," + Double.toString(value.get(hapIndex, readIndex)));
            }
            pw.println("");
        }
        pw.close();

        // add entry
        addEntry(loc, name, os.toByteArray());
    }

    public void addReads(final Locatable loc, final String name, final List<GATKRead> value, List<GATKRead> readsList1, List<GATKRead> readsList2, final String msg) throws IOException {

        // build text representation
        final ByteArrayOutputStream   os = new ByteArrayOutputStream();
        final PrintWriter             pw = new PrintWriter(os);
        if (WIDE_READS) {
            final List<String> names1 = readsList1.stream().map(GATKRead::getName).collect(Collectors.toList());
            final List<String> names2 = readsList2.stream().map(GATKRead::getName).collect(Collectors.toList());
            pw.println("name,supp,l1,l2,msg");
            for (GATKRead read : value) {
                final String readName = read.getName();
                pw.println(String.format("%s,%d,%d,%d,%s",
                        readName,
                        read.isSupplementaryAlignment() ? 1 : 0,
                        names1.contains(readName) ? 1 : 0,
                        names2.contains(readName) ? 1 : 0,
                        msg
                        )
                );
            }
        } else {
            pw.println("name,supp");
            for (GATKRead read : value) {
                pw.println(read.getName() + "," + (read.isSupplementaryAlignment() ? 1 : 0));
            }
        }
        pw.close();

        // add entry
        addEntry(loc, name, os.toByteArray());
    }
}
