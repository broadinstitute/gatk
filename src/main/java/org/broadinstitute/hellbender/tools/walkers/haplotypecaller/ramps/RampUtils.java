package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Collection;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

public class RampUtils {

    protected static final Logger logger = LogManager.getLogger(RampUtils.class);

    public static class GATKReadComparator implements Comparator<GATKRead> {

        @Override
        public int compare(final GATKRead o1, final GATKRead o2) {

            int     delta = (o1.isReverseStrand() ? 1 : 0) - (o2.isReverseStrand() ? 1 : 0);

            if ( delta == 0 ) {
                delta = o1.commonToString().compareTo(o2.commonToString());
                if ( delta == 0 ) {
                    delta = o1.getBasesString().compareTo(o2.getBasesString());
                    if ( delta == 0 ) {
                        delta = (new String(o1.getBaseQualitiesNoCopy())).compareTo((new String(o2.getBaseQualitiesNoCopy())));
                        if ( delta == 0 ) {
                            delta = o1.getSoftStart() - o2.getSoftStart();
                            if ( delta == 0 ) {
                                delta = o1.getSoftEnd() - o2.getSoftEnd();
                                if ( delta == 0 ) {
                                    delta = o1.getStart() - o2.getStart();
                                    if ( delta == 0 ) {
                                        delta = o1.getEnd() - o2.getEnd();
                                        if ( delta == 0 ) {
                                            delta = o1.getUnclippedStart() - o2.getUnclippedStart();
                                            if ( delta == 0 ) {
                                                delta = o1.getUnclippedEnd() - o2.getUnclippedEnd();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return delta;
        }
    }

    public static class HaplotypeComparator implements Comparator<Haplotype> {

        @Override
        public int compare(final Haplotype o1, final Haplotype o2) {

            int         delta = (o1.isReference() ? 1 : 0) - (o2.isReference() ? 1 : 0);

            if ( delta == 0 ) {
                delta = o1.getLocation().getContig().compareTo(o2.getLocation().getContig());
                if ( delta == 0 ) {
                    delta = o1.getLocation().getStart() - o2.getLocation().getStart();
                    if ( delta == 0 ) {
                        delta = o1.getLocation().getEnd() - o2.getLocation().getEnd();
                        if ( delta == 0 ) {
                            double      ddelta = o1.getScore() - o2.getScore();
                            if ( ddelta < 0 )
                                delta = -1;
                            else if ( ddelta > 0 )
                                delta = 1;
                            else
                                delta = 0;
                            if ( delta == 0 ) {
                                delta = o1.getBaseString().compareTo(o2.getBaseString());
                            }
                        }
                    }
                }
            }

            return delta;
        }
    }

    static public void compareHaplotypes(final Collection<Haplotype> c1, final Collection<Haplotype> c2) {

        // must be the same length
        if ( c1.size() != c2.size() )
            throw new RuntimeException("haplotype size verification failed: " +  c1.size() + " vs " + c2.size());

        final LinkedList<Haplotype>     l1 = new LinkedList<>(c1);
        final LinkedList<Haplotype>     l2 = new LinkedList<>(c2);

        // sort and compare
        final HaplotypeComparator       hc = new HaplotypeComparator();
        l1.sort(hc);
        l2.sort(hc);
        for ( int i = 0 ; i < l1.size() ; i++ ) {
            if ( hc.compare(l1.get(i), l2.get(i)) != 0 ) {
                throw new RuntimeException("haplotype failed verification on index " + i);
            }
        }
    }

    public static void compareReads(final Collection<GATKRead> c1, final Collection<GATKRead> c2) {

        // must be the same length
        if ( c1.size() != c2.size() )
            throw new RuntimeException("reads size verification failed: " +  c1.size() + " vs " + c2.size());

        final LinkedList<GATKRead>       l1 = new LinkedList<>(c1);
        final LinkedList<GATKRead>       l2 = new LinkedList<>(c2);

        // sort and compare
        final GATKReadComparator         rc = new GATKReadComparator();
        l1.sort(rc);
        l2.sort(rc);
        for ( int i = 0 ; i < l1.size() ; i++ ) {
            if ( rc.compare(l1.get(i), l2.get(i)) != 0 ) {
                rc.compare(l1.get(i), l2.get(i));
                logger.error("l1 read: " + readInfo(l1.get(i)));
                logger.error("l2 read: " + readInfo(l2.get(i)));
                throw new RuntimeException("reads failed verification on index " + i);
            }
        }
    }

    public static String readInfo(final GATKRead read) {
        final String      prefix = read.isSupplementaryAlignment() ? "1," : "0,";
        if ( read instanceof SAMRecordToGATKReadAdapter) {
            return prefix + ((SAMRecordToGATKReadAdapter) read).getEncapsulatedSamRecord().toString();
        } else {
            return prefix + read.toString();
        }
    }

    public static void logReads(final String debugReadsParam, final String message, final Collection<GATKRead> reads) {
        if ( debugReadsParam == null ) {
            return;
        }

        // filter?
        String debugReads = debugReadsParam;
        if ( debugReads.contains("|") ) {
            String[]    toks = debugReads.split("\\|");
            if ( toks.length > 1 ) {
                final String      filter = toks[0];
                debugReads = toks[1];

                if ( !message.contains(filter) )
                    return;
            }
        }

        logger.info(message + ", " + reads.size() + " reads");

        List<String>        readNames = null;
        if ( !debugReads.equals("ALL") ) {
            readNames = new LinkedList<>();
            for ( String name : debugReads.split(",") )
                readNames.add(name);
        }

        for ( GATKRead read : reads ) {
            if ( readNames == null || readNames.contains(read.getName()) )
                logger.info(readInfo(read));
        }
    }
}
