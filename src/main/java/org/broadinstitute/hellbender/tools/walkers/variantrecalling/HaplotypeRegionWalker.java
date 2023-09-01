package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * a service class for HaplotypeBasedVariableRecaller that reads a SAM/BAM file, interprets the reads as haplotypes
 * and called a provided consumer with the 'best' haplotypes found for a given query location.
 *
 * The notion of best is controlled by a fitnessScore - which is related to the distance between the query location
 * and the haplotype boundaries (less is better)
 */
public class HaplotypeRegionWalker {

    private static final Logger logger = LogManager.getLogger(HaplotypeRegionWalker.class);

    final private SamReader             samReader;
    final private List<Haplotype> walkerHaplotypes = new LinkedList<>();
    final private boolean         reusePreviousResults = false;

    HaplotypeRegionWalker(final HaplotypeBasedVariantRecallerArgumentCollection vrArgs, final Path referencePath, final int cloudPrefetchBuffer) {
        final Path samPath = IOUtils.getPath(vrArgs.haplotypesBamFile);

        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
        samReader = SamReaderFactory.makeDefault().referenceSequence(referencePath).open(samPath, cloudWrapper, cloudIndexWrapper);
    }

    // iterate over a haplotypes behind the variant context. If there are multiple sets, select the best set
    // see fitnessScore for the scoring function
    void forBest(final Locatable queryLoc, final Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);

        final List<Haplotype> best = new LinkedList<>();

        forEach(queryLoc, haplotypes -> {
            if ( best.size() == 0 ||
                    (fitnessScore(queryLoc, haplotypes) > fitnessScore(queryLoc, best)) ) {
                best.clear();
                best.addAll(haplotypes);
            }
        });

        if ( best.size() != 0 ) {
            action.accept(best);
        }
    }

    private double fitnessScore(final Locatable loc, final List<Haplotype> haplotypes) {
        Objects.requireNonNull(haplotypes);
        if ( haplotypes.size() == 0 ) {
            return 0;
        }
        final Locatable   hloc = haplotypes.get(0).getGenomeLocation();

        // determine spacing before and end of loc
        final int         before = Math.max(1, loc.getStart() - hloc.getStart());
        final int         after = Math.max(1, hloc.getEnd() - loc.getEnd());

        // score reflects closeness to being in the center
        final double       score = 1.0 - 2 * Math.abs(0.5 - (double)before / (before + after));

        if ( logger.isDebugEnabled() ) {
            logger.debug(String.format("loc %s, hloc: %s, before: %d, after: %d, score: %f",
                    loc, hloc, before, after, score));
        }

        return score;
    }


    void forEach(final Locatable queryLoc, final Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);

        // use last results? (note that this can be problematic if a vc is inside two haplotype areas)
        if ( reusePreviousResults
                && (walkerHaplotypes.size() != 0)
                && walkerHaplotypes.get(0).getGenomeLocation().contains(queryLoc) ) {
            action.accept(walkerHaplotypes);
        } else {

            // must query
            walkerHaplotypes.clear();
            final SAMRecordIterator iter = samReader.query(queryLoc.getContig(), queryLoc.getStart(), queryLoc.getEnd(), false);
            iter.forEachRemaining(record -> {
                if (isHaplotypeRecord(record)) {
                    Locatable loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                    if ((walkerHaplotypes.size() != 0) && !walkerHaplotypes.get(0).getGenomeLocation().equals(loc)) {
                        action.accept(walkerHaplotypes);
                        walkerHaplotypes.clear();
                    }
                    walkerHaplotypes.add(buildHaplotype(record));
                }
            });
            if (walkerHaplotypes.size() != 0) {
                action.accept(walkerHaplotypes);
            }
            iter.close();
        }
    }

    private boolean isHaplotypeRecord(final SAMRecord record) {
        return record.getReadName().startsWith("HC_");
    }

    private Haplotype buildHaplotype(final SAMRecord record) {

        final Locatable loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
        final Haplotype haplotype = new Haplotype(record.getReadBases(), loc);
        haplotype.setCigar(record.getCigar());
        return haplotype;
    }
}
