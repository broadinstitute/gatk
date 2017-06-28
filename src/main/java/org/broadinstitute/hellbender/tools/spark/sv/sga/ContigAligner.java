package org.broadinstitute.hellbender.tools.spark.sv.sga;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

public class ContigAligner {

    private final String indexImageFile;

    ContigAligner(final String indexImageFile) {
        this.indexImageFile = indexImageFile;
    }

    /**
     * Takes a collection of assembled contigs and aligns them to the reference with bwa-mem. Non-canonical
     * (secondary) alignments are filtered out, preserving the primary and supplementary alignments.
     * Within the output list, alignments are sorted first by contig (based on the order in which
     * the contigs were passed in, and then by their start position on the contig).
     * @param assemblyId An identifier for the assembly or set of contigs
     * @param contigsCollection The set of all canonical (primary or supplementary) alignments for the contigs.*/

    AlignedAssembly alignContigs(final int assemblyId, final ContigsCollection contigsCollection) {

        final List<AlignedContig> alignedContigs = new ArrayList<>(contigsCollection.getContents().size());

        final BwaMemIndex index = BwaMemIndexSingleton.getInstance(indexImageFile);

        try ( final BwaMemAligner aligner = new BwaMemAligner(index) ) {

            final List<String> refNames = index.getReferenceContigNames();

            final List<Tuple2<ContigsCollection.ContigID, ContigsCollection.ContigSequence>> contents = contigsCollection.getContents();

            final List<byte[]> seqs = contents.stream().map(contigInfo -> contigInfo._2.toString().getBytes()).collect(Collectors.toList());

            final List<List<BwaMemAlignment>> allAlignments = aligner.alignSeqs(seqs);
            for (int contigIdx = 0; contigIdx < seqs.size(); ++contigIdx) {

                final int contigLen = seqs.get(contigIdx).length;

                // filter out secondary alignments, convert to AlignmentInterval objects and sort by alignment start pos
                final List<AlignedAssembly.AlignmentInterval> alignmentIntervals
                        = allAlignments.get(contigIdx).stream()
                        .filter(a -> (a.getSamFlag()&SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue())==0)
                        .filter(a -> (a.getSamFlag()&SAMFlag.READ_UNMAPPED.intValue())==0)
                        .map(a -> new AlignedAssembly.AlignmentInterval(a, refNames, contigLen))
                        .sorted(Comparator.comparing(a -> a.startInAssembledContig))
                        .collect(Collectors.toList());

                final String contigName =
                        AlignedAssemblyOrExcuse.formatContigName(assemblyId,
                                Integer.valueOf(contents.get(contigIdx)._1.toString().split(" ")[0].replace("contig", "").replace("-", "").replace(">", "")));

                alignedContigs.add( new AlignedContig(contigName, seqs.get(contigIdx), alignmentIntervals) );
            }
        }

        return new AlignedAssembly(assemblyId, alignedContigs);
    }
}
