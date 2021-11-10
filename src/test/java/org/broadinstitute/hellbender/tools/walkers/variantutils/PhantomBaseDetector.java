package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.table.TableFeature;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Walk vcf and look at phase shifts in hmers.",
        summary = "Walk vcf and look at phase shifts in hmers.",
        programGroup = ReferenceProgramGroup.class
)
public class PhantomBaseDetector extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private GATKPath outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF or BED file)")
    public FeatureInput<TableFeature> featuresFile = null;

    private PrintStream outputStream = null;


    @Override
    public void onTraversalStart() {
        outputStream = outputFile != null ? new PrintStream(outputFile.getOutputStream()) : System.out;
        outputStream.println("HEADER GT GQ MaxHmerLength Call");
    }


    /** Variables used for counting reads of various types */
    private int snps_on_hmers_edge = 0;
    private int insertions_on_hmers_edge = 0;
    private int deletions_on_hmers = 0;
    private int non_vars_on_hmers = 0;

    // SNP counting variables; first column for phantom bases, second others
    private int[] AtoC = new int[] {0,0};
    private int[] AtoG = new int[] {0,0};
    private int[] AtoT = new int[] {0,0};
    private int[] CtoA = new int[] {0,0};
    private int[] CtoG = new int[] {0,0};
    private int[] CtoT = new int[] {0,0};
    private int[] GtoA = new int[] {0,0};
    private int[] GtoC = new int[] {0,0};
    private int[] GtoT = new int[] {0,0};
    private int[] TtoA = new int[] {0,0};
    private int[] TtoC = new int[] {0,0};
    private int[] TtoG = new int[] {0,0};

    private void calcSNP(char a, char b, boolean isPhantom) {
        switch(a) {
            case 'A':
                switch(b) {
                    case 'C':
                        if (isPhantom) {
                            AtoC[0]++;
                        } else {
                            AtoC[1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            AtoG[0]++;
                        } else {
                            AtoG[1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            AtoT[0]++;
                        } else {
                            AtoT[1]++;
                        }
                        break;
                }
                break;
            case 'C':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            CtoA[0]++;
                        } else {
                            CtoA[1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            CtoG[0]++;
                        } else {
                            CtoG[1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            CtoT[0]++;
                        } else {
                            CtoT[1]++;
                        }
                        break;
                }
                break;
            case 'G':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            GtoA[0]++;
                        } else {
                            GtoA[1]++;
                        }
                        break;
                    case 'C':
                        if (isPhantom) {
                            GtoC[0]++;
                        } else {
                            GtoC[1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            GtoT[0]++;
                        } else {
                            GtoT[1]++;
                        }
                        break;
                }
                break;
            case 'T':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            TtoA[0]++;
                        } else {
                            TtoA[1]++;
                        }
                        break;
                    case 'C':
                        if (isPhantom) {
                            TtoC[0]++;
                        } else {
                            TtoC[1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            TtoG[0]++;
                        } else {
                            TtoG[1]++;
                        }
                        break;
                }
                break;
            default:
        }
    }

    // Takes in relative ref position and cigar list to get relative read position when walking
    private int getRelativeReadPosition(int relativeRefPosition, List<CigarOperator> cigarList) {
        int relativeReadPosition = 0;
        int insertCount = 0;
        for (int i = 0; i < relativeRefPosition; i++) {
            CigarOperator c = cigarList.get(i);
            switch (c) {
                case D:
                    break;
                case I:
                    insertCount++;
                    break;
                default:
                    relativeReadPosition += 1 + insertCount;
                    insertCount = 0;
            }
        }
        return relativeReadPosition;
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        Cigar cigar = read.getCigar();

        // Get cigar as list of operators
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        List<CigarOperator> cigarList = new ArrayList<>();
        for (CigarElement currentCigar : cigarElements) {
            for (int j = 0; j < currentCigar.getLength(); j++) {
                cigarList.add(currentCigar.getOperator());
            }
        }

        // Loop over all hmer intervals overlapping current read
        if ( featureContext.getValues(featuresFile).isEmpty() ) {
            String test = "hi";
        }
        if ( ! featureContext.getValues(featuresFile).isEmpty() ) {
            for (TableFeature hmer : featureContext.getValues(featuresFile)) {
                // Get interval for hmer from table
                final SimpleInterval interval = new SimpleInterval(hmer.get("HEADER"));

                // (Relative, e.g. 0 = start of read) Start & End location of hmer interval intersected with read interval
                final int relRefIntStart = Math.max(interval.getStart() - read.getStart(), 0);
                final int relRefIntEnd = relRefIntStart + Math.min(interval.getEnd() - interval.getStart(), read.getEnd() - interval.getStart());
                final int readOverlapLength = relRefIntEnd - relRefIntStart;

                // Get relative position read should start on clipping edge overlap cases
                final int relReadIntStart = getRelativeReadPosition(relRefIntStart, cigarList);
                final int relReadIntEnd = getRelativeReadPosition(relRefIntEnd, cigarList);

                // Shorten list of bases and cigar list to just focus on overlap
                final byte[] readBasesOverlap = Arrays.copyOfRange(read.getBases(), relReadIntStart, relReadIntEnd);
                final List<CigarOperator> cigarOverlap = cigarList.subList(relReadIntStart, relReadIntEnd);

                final int hmerLength = Integer.parseInt(hmer.get("length"));
                final char hmerBase = hmer.get("base").charAt(0);

                final boolean edgeOverlap = ( read.getEnd() < interval.getEnd() ) || ( read.getStart() > interval.getStart() );

                for (int i = 0; i < readOverlapLength; i++) {
                    final char readBase = (char) readBasesOverlap[i];
                    if (readBase != hmerBase) {
                        switch (cigarOverlap.get(i)) {
                            case S:
                                break;

                            case M:
                                if (edgeOverlap) {
                                    snps_on_hmers_edge++;
                                } else {
                                    char edgeBase;
                                    // Determine base at end of hmer based on strand
                                    if (read.isReverseStrand()) {
                                        edgeBase = (char) read.getBase(relReadIntStart - 1);
                                    } else {
                                        edgeBase = (char) read.getBase(relReadIntEnd + 1);
                                    }

                                    // Check if edge base agrees with SNP base
                                    boolean phantomBase;
                                    if (edgeBase == readBase) {
                                        phantomBase = true;
                                    } else {
                                        phantomBase = false;
                                    }

                                    // Increment appropriate counter for SNP
                                    calcSNP(hmerBase, readBase, phantomBase);
                                }
                                break;

                            default:
                                break;
                        }
                    }
                }

                String test = "Hello World";
            }
        }
    }


    @Override
    public void closeTool() {
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
