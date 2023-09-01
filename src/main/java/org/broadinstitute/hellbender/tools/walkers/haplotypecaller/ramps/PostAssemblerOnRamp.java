package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LongHomopolymerHaplotypeCollapsingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;
import java.util.Optional;

public class PostAssemblerOnRamp extends OnRampBase {

    final private Path            haplotypeBAMWriterPath;
    final private Path            haplotypeBAIWriterPath;

    final private SamReader       haplotypeReader;

    public PostAssemblerOnRamp(final String filename) throws IOException {
        super(filename);

        // open haplotype file (bam)
        haplotypeBAMWriterPath = File.createTempFile("haplotypes_", ".bam").toPath();
        haplotypeBAIWriterPath = getBamIndexPath(haplotypeBAMWriterPath);
        copyStreamToPath(getEntry(null, "haplotypes.bam"), haplotypeBAMWriterPath);
        copyStreamToPath(getEntry(null, "haplotypes.bai"), haplotypeBAIWriterPath);
        final SamInputResource    samInputResource = SamInputResource.of(haplotypeBAMWriterPath);
        samInputResource.index(haplotypeBAIWriterPath);
        haplotypeReader = SamReaderFactory.makeDefault().open(samInputResource);
    }

    @Override
    public void close() throws IOException {
        if ( haplotypeReader != null ) {
            haplotypeReader.close();
            haplotypeBAMWriterPath.toFile().delete();
            haplotypeBAIWriterPath.toFile().delete();
        }
        super.close();
    }

    public Optional<AssemblyRegion> getOptionalAssemblyRegion(final AssemblyRegion region, final String name, final SAMFileHeader header) throws IOException {

        final JSONObject                  obj = getRegion(region.getSpan(), name);
        if ( obj == null ) {
            return Optional.empty();
        } else {
            return Optional.of(parseAssemblyRegion(obj, header));
        }
    }

    private AssemblyRegion parseAssemblyRegion(final JSONObject obj, final SAMFileHeader header) {
        return new AssemblyRegion(
                new SimpleInterval(obj.getString("activeSpan")),
                new SimpleInterval(obj.getString("paddedSpan")),
                true,
                header);
    }

    public AssemblyResultSet getAssembyResult(final AssemblyRegion region, final String name, final SAMFileHeader header, final HaplotypeCallerArgumentCollection hcArgs,
                                              final Logger logger, final SmithWatermanAligner aligner, final ReferenceSequenceFile referenceReader) throws IOException {

        // create result
        final AssemblyResultSet       assemblyResult = new AssemblyResultSet();

        // read fields
        final JSONObject              obj = getRegion(region.getSpan(), name);
        if ( obj == null )
            return null;
        final List<Haplotype>         haplotypes = readHaplotypes(region);
        final SimpleInterval          paddedReferenceLoc = new SimpleInterval(obj.getString("paddedReferenceLoc"));
        final AssemblyRegion          regionForGenotyping = parseAssemblyRegion(obj.getJSONObject("regionForGenotyping"), header);

        // DEBUG, verify haplotypes
        /*
        boolean                 debugVerifyHaplotypes = true;
        if ( debugVerifyHaplotypes ) {
            List<Haplotype>     verifyHaplotypes = getHaplotypes(region, name + ".haplotypes");

            RampUtils.compareHaplotypes(haplotypes, verifyHaplotypes);
        }
        */

        // install
        assemblyResult.setPaddedReferenceLoc(paddedReferenceLoc);
        assemblyResult.setRegionForGenotyping(regionForGenotyping);
        haplotypes.forEach(h -> assemblyResult.add(h));

        // assemblyResult.getFullReferenceWithPadding()
        final ReferenceSequence       refSeq = referenceReader.getSubsequenceAt(paddedReferenceLoc.getContig(),
                                                paddedReferenceLoc.getStart(), paddedReferenceLoc.getEnd());
        assemblyResult.setFullReferenceWithPadding(refSeq.getBases());

        // refView
        final LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing = (hcArgs.flowAssemblyCollapseHKerSize > 0 && LongHomopolymerHaplotypeCollapsingEngine.needsCollapsing(assemblyResult.getReferenceHaplotype().getBases(), hcArgs.flowAssemblyCollapseHKerSize, logger))
                ? new LongHomopolymerHaplotypeCollapsingEngine(hcArgs.flowAssemblyCollapseHKerSize, hcArgs.flowAssemblyCollapsePartialMode, assemblyResult.getFullReferenceWithPadding(),
                paddedReferenceLoc, logger, hcArgs.assemblerArgs.debugAssembly, aligner, hcArgs.getHaplotypeToReferenceSWParameters())
                : null;
        assemblyResult.setHaplotypeCollapsingEngine(haplotypeCollapsing);

        return assemblyResult;
    }

    private List<Haplotype> readHaplotypes(final AssemblyRegion region) {

        // read records within the region range
        final List<Haplotype>     haplotypes = new LinkedList<>();
        final Locatable           loc = region;
        final String              locStr = String.format("%s:%d-%d", loc.getContig(), loc.getStart(), loc.getEnd());
        final SAMRecordIterator   iter = haplotypeReader.query(loc.getContig(), loc.getStart(), loc.getEnd(), false);
        while ( iter.hasNext() ) {
            final SAMRecord   record = iter.next();

            // get special attribute
            boolean         isReference = false;
            double          score = Double.NaN;
            int             alignmentStartHapwrtRef = 0;
            final String          specialAttr = record.getStringAttribute(AssemblyBasedCallerUtils.EXT_SPECIAL_TAG);
            if ( specialAttr != null ) {
                String[]      specialAttrToks = specialAttr.split(",");
                isReference = (Integer.parseInt(specialAttrToks[0]) != 0);
                score = Double.parseDouble(specialAttrToks[1]);
                alignmentStartHapwrtRef = Integer.parseInt(specialAttrToks[2]);

                // filter?
                if ( !specialAttrToks[3].equals(locStr) )
                    continue;;
            }
            final Haplotype       haplotype = new Haplotype(record.getReadBases(), isReference);
            haplotype.setGenomeLocation(new SimpleInterval(record.getContig(), record.getStart(), record.getEnd()));
            haplotype.setCigar(record.getCigar());
            haplotype.setAlignmentStartHapwrtRef(alignmentStartHapwrtRef);
            if ( !Double.isNaN(score) )
                haplotype.setScore(score);

            haplotypes.add(haplotype);
        }
        iter.close();

        return haplotypes;
    }
}
