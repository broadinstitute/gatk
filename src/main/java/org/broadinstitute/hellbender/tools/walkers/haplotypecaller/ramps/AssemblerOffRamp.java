package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.json.JSONObject;

import java.io.*;
import java.nio.file.Path;
import java.util.*;

public class AssemblerOffRamp extends OffRampBase {

    private Path                haplotypeBAMWriterPath;
    private HaplotypeBAMWriter  haplotypeBAMWriter;

    public AssemblerOffRamp(String filename) throws IOException {
        super(filename);
    }

    @Override
    public void close() throws IOException {
        if ( haplotypeBAMWriter != null ) {
            haplotypeBAMWriter.close();
            addEntry(null, "haplotypes.bam", haplotypeBAMWriterPath);
            addEntry(null, "haplotypes.bai", getBamIndexPath(haplotypeBAMWriterPath));
            haplotypeBAMWriterPath.toFile().delete();
        }
        super.close();
    }

    public synchronized void add(final Locatable loc, final String name, final AssemblyResultSet value,
                                 final Optional<AssemblyRegion> nonVariantLeftFlankRegion,
                                 final Optional<AssemblyRegion> nonVariantRightFlankRegion,
                                 final SAMFileHeader header) throws IOException {

        // add region
        final JSONObject regionObj = new JSONObject();
        regionObj.put("contig", loc.getContig());
        regionObj.put("start", loc.getStart());
        regionObj.put("end", loc.getEnd());
        final JSONObject assemblerObj = new JSONObject();

        if ( value != null ) {
            assemblerObj.put("paddedReferenceLoc", value.getPaddedReferenceLoc());
            assemblerObj.put("regionForGenotyping", getJson(value.getRegionForGenotyping()));
        }
        if ( nonVariantLeftFlankRegion != null ) {
            nonVariantLeftFlankRegion.ifPresent(v -> assemblerObj.put("nonVariantLeftFlankRegion", getJson(v)));
        }
        if ( nonVariantRightFlankRegion != null ) {
            nonVariantRightFlankRegion.ifPresent(v -> assemblerObj.put("nonVariantRightFlankRegion", getJson(v)));
        }

        // add haplotypes
        if ( value != null ) {
            //addHaplotypes(loc, name + ".haplotypes", value.getHaplotypeList());
            writeHaplotypes(loc, value, header);
        }

        // add to regions
        if ( assemblerObj.length() != 0 ) {
            regionObj.put(name, assemblerObj);
        }
        info.getJSONArray("regions").put(regionObj);
    }

    private synchronized void writeHaplotypes(final Locatable loc, final AssemblyResultSet assemblyResultSet, final SAMFileHeader header) throws IOException {

        // TEMP logic
        // open haplotype writer?
        if ( haplotypeBAMWriter == null ) {

            haplotypeBAMWriterPath = File.createTempFile("haplotypes_", ".bam").toPath();
            haplotypeBAMWriter = new HaplotypeBAMWriter(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES,
                                            haplotypeBAMWriterPath,
                                            true,
                                            false,
                                            header
                                            );
            haplotypeBAMWriter.setAddSpecialTag(true);
        }

        haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                assemblyResultSet.getHaplotypeList(),
                assemblyResultSet.getPaddedReferenceLoc(),
                Collections.emptySet(),
                Collections.emptySet(),
                null,
                loc);
    }

    private JSONObject getJson(final AssemblyRegion value) {

        final JSONObject      json = new JSONObject();

        json.put("activeSpan", value.getSpan());
        json.put("paddedSpan", value.getPaddedSpan());

        return json;
    }
}