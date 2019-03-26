package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.sam.util.PhysicalLocationInt;

import picard.sam.markduplicates.MarkDuplicates;

public class ReadElementExtractorImpl {
    static class Tile implements ReadElementExtractor {
        @Override public String header() { return "TILE"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractTile(read));
        }
    }

    static class XCoord implements ReadElementExtractor {
        @Override public String header() { return "X_COORD"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractX(read));
        }
    }

    static class YCoord implements ReadElementExtractor {
        @Override public String header() { return "Y_COORD"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractY(read));
        }
    }

    static class InsertSize implements ReadElementExtractor {
        @Override public String header() { return "I_SIZE"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(read.getFragmentLength());
        }
    }

    static class Duplicate implements ReadElementExtractor {
        @Override public String header() { return "DUP"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return extractDuplicateType(read);
        }
    }

    static class Errors implements ReadElementExtractor {
        @Override public String header() { return "NM"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(read.getAttributeAsInteger("NM"));
        }
    }

    static class Length implements ReadElementExtractor {
        @Override public String header() { return "ALIGNED_LENGTH"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractNonclippedLength(read));
        }
    }

    static class Length2 implements ReadElementExtractor {
        @Override public String header() { return "HIGH_QUAL_LENGTH"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractNonLowQualLength(read));
        }
    }


    static class MappingQ implements ReadElementExtractor {
        @Override public String header() { return "MAPPING_Q"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(read.getMappingQuality());
        }
    }

    static class Mapped implements ReadElementExtractor {
        @Override public String header() { return "MAPPED"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(!read.isUnmapped());
        }
    }
    static class First implements ReadElementExtractor {
        @Override public String header() { return "FIRST"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(read.isFirstOfPair());
        }
    }

    static class BaseQual implements ReadElementExtractor {
        @Override public String header() { return "BASE_QUAL"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return String.valueOf(extractMeanQual(read));
        }
    }

    static class ReadGroup implements ReadElementExtractor {
        @Override public String header() { return "READGROUP"; }
        @Override public String extractElement(final GATKRead read, final SAMFileHeader header) {
            return ReadUtils.getPlatformUnit(read, header);
        }
    }


    //average base quality

    private static OpticalDuplicateFinder opticalDuplicateFinder = new OpticalDuplicateFinder();

    private static short extractTile(final GATKRead read) {
        final PhysicalLocationInt location = new PhysicalLocationInt();
        opticalDuplicateFinder.addLocationInformation(read.getName(),location);
        return location.tile;
    }

    private static int extractX(final GATKRead read) {
        final PhysicalLocationInt location = new PhysicalLocationInt();
        opticalDuplicateFinder.addLocationInformation(read.getName(),location);
        return location.x;
    }

    private static int extractY(final GATKRead read) {
        final PhysicalLocationInt location = new PhysicalLocationInt();
        opticalDuplicateFinder.addLocationInformation(read.getName(),location);
        return location.y;
    }

   private static String extractDuplicateType(final GATKRead read){
        if (!read.isDuplicate()) return "NOT";
        if (!read.hasAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG)) return "DUP";

        return read.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG);
   }

   private static int extractNonclippedLength(final GATKRead read){
       if (read.isUnmapped()) return 0;
        return read.getCigar().getCigarElements().stream()
                .filter(e->e.getOperator()!= CigarOperator.HARD_CLIP && e.getOperator()!=CigarOperator.SOFT_CLIP)
                .filter(e->e.getOperator().consumesReadBases())
                .mapToInt(CigarElement::getLength)
                .sum();
   }

    private static int extractNonLowQualLength(final GATKRead read) {
        int count = 0;
        for (int i = 0; i < read.getBaseQualityCount(); i++) {
            if (read.getBaseQuality(i) > 2) {
                count++;
            }
        }
        return count;
    }

    private static double extractMeanQual(final GATKRead read) {
        int count = 0;
        for (int i = 0; i < read.getBaseQualityCount(); i++) {
                count++;
        }
        return count/(double)read.getBaseQualities().length;
    }
}
