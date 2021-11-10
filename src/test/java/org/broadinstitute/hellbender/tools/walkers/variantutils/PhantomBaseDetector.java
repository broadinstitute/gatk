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
        oneLineSummary = "Walk vcf and look at phantom bases in hmers.",
        summary = "Walk vcf and look at phantom bases in hmers.",
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
//        outputStream.println("HEADER GT GQ MaxHmerLength Call");
    }


    /** Variables used for counting reads of various types */
    private int snps_on_hmers_edge = 0;
    private int insertions_on_hmers_edge = 0;
    private int deletions_on_hmers = 0;
    private int non_vars_on_hmers = 0;

    // SNP counting variables; first column for phantom bases, second others
    private int[][] AtoC = new int[32][2];
    private int[][] AtoG = new int[32][2];
    private int[][] AtoT = new int[32][2];
    private int[][] CtoA = new int[32][2];
    private int[][] CtoG = new int[32][2];
    private int[][] CtoT = new int[32][2];
    private int[][] GtoA = new int[32][2];
    private int[][] GtoC = new int[32][2];
    private int[][] GtoT = new int[32][2];
    private int[][] TtoA = new int[32][2];
    private int[][] TtoC = new int[32][2];
    private int[][] TtoG = new int[32][2];

    private void calcSNP(char a, char b, boolean isPhantom, int hmer_position) {
        // Group together all changes happening past 30 bases into hmer
        hmer_position = Math.min(hmer_position, 31);
        switch(a) {
            case 'A':
                switch(b) {
                    case 'C':
                        if (isPhantom) {
                            AtoC[hmer_position][0]++;
                        } else {
                            AtoC[hmer_position][1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            AtoG[hmer_position][0]++;
                        } else {
                            AtoG[hmer_position][1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            AtoT[hmer_position][0]++;
                        } else {
                            AtoT[hmer_position][1]++;
                        }
                        break;
                }
                break;
            case 'C':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            CtoA[hmer_position][0]++;
                        } else {
                            CtoA[hmer_position][1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            CtoG[hmer_position][0]++;
                        } else {
                            CtoG[hmer_position][1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            CtoT[hmer_position][0]++;
                        } else {
                            CtoT[hmer_position][1]++;
                        }
                        break;
                }
                break;
            case 'G':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            GtoA[hmer_position][0]++;
                        } else {
                            GtoA[hmer_position][1]++;
                        }
                        break;
                    case 'C':
                        if (isPhantom) {
                            GtoC[hmer_position][0]++;
                        } else {
                            GtoC[hmer_position][1]++;
                        }
                        break;
                    case 'T':
                        if (isPhantom) {
                            GtoT[hmer_position][0]++;
                        } else {
                            GtoT[hmer_position][1]++;
                        }
                        break;
                }
                break;
            case 'T':
                switch(b) {
                    case 'A':
                        if (isPhantom) {
                            TtoA[hmer_position][0]++;
                        } else {
                            TtoA[hmer_position][1]++;
                        }
                        break;
                    case 'C':
                        if (isPhantom) {
                            TtoC[hmer_position][0]++;
                        } else {
                            TtoC[hmer_position][1]++;
                        }
                        break;
                    case 'G':
                        if (isPhantom) {
                            TtoG[hmer_position][0]++;
                        } else {
                            TtoG[hmer_position][1]++;
                        }
                        break;
                }
                break;
            default:
                break;
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

//        if ( featureContext.getValues(featuresFile).isEmpty() ) {
//            String test = "hi";
//        }
//        if ( ! featureContext.getValues(featuresFile).isEmpty() ) {
        for (TableFeature hmer : featureContext.getValues(featuresFile)) {
            // Get interval for hmer from reference table
            final SimpleInterval hmer_interval = new SimpleInterval(hmer.get("HEADER"));

            // (Relative, e.g. 0 = start of read) Start & End location of hmer interval intersected with read interval
            final int relRefIntStart = Math.max(hmer_interval.getStart() - read.getStart(), 0);
            final int relRefIntEnd = relRefIntStart +
                    Math.min(hmer_interval.getEnd() - hmer_interval.getStart(),
                            Math.min(read.getEnd() - hmer_interval.getStart(), read.getEnd() - read.getStart()));

            // Get relative position read should start on clipping edge overlap cases
            final int relReadIntStart = getRelativeReadPosition(relRefIntStart, cigarList);
            final int relReadIntEnd = getRelativeReadPosition(relRefIntEnd, cigarList);
            final int readOverlapLength = relReadIntEnd - relReadIntStart;

            // Shorten list of bases and cigar list to just focus on overlap
            final byte[] readBasesOverlap = Arrays.copyOfRange(read.getBases(), relReadIntStart, relReadIntEnd);
            final List<CigarOperator> cigarOverlap = cigarList.subList(relReadIntStart, relReadIntEnd);

            final int hmerLength = Integer.parseInt(hmer.get("length"));
            final char hmerBase = hmer.get("base").charAt(0);

            final boolean edgeOverlap = ( read.getEnd() <= hmer_interval.getEnd() ) || ( read.getStart() >= hmer_interval.getStart() );

            // Loop over intersected read and compare bases against hmer
            for (int i = 0; i < readOverlapLength; i++) {
                final char readBase = (char) readBasesOverlap[i];
                boolean seen_snp = false;
                if (readBase != hmerBase) {
                    switch (cigarOverlap.get(i)) {
                        case S:
                            break;

                        case M:
                            if (edgeOverlap) {
                                snps_on_hmers_edge++;
                            } else {
                                char edgeBase;
                                int relBasePosition; // in strand direction, 0 indexed
                                final int unseen_hmers = Math.max(hmerLength - readOverlapLength, 0); // num hmers not seen in overlap from left

                                // Determine base at end of hmer based on strand
                                if (read.isReverseStrand()) {
                                    edgeBase = (char) referenceContext.getBases()[hmer.getStart()-referenceContext.getStart()-1];
                                    relBasePosition = unseen_hmers + readOverlapLength - i;
                                } else {
                                    edgeBase = (char) referenceContext.getBases()[hmer.getEnd()-referenceContext.getStart()+1];
                                    //edgeBase = (char) read.getBase(relReadIntEnd + 1);
                                    relBasePosition = unseen_hmers + i;
                                }

                                // Check if edge base agrees with SNP base
                                boolean phantomBase;
                                if (edgeBase == readBase) {
                                    phantomBase = true;
                                } else {
                                    phantomBase = false;
                                }

                                // Increment appropriate counter for SNP
                                calcSNP(hmerBase, readBase, phantomBase, relBasePosition);
                                seen_snp = true;
                            }
                            break;

                        case I:
                            break;

                        case D:
                            break;

                        default:
                            break;
                    }
                } else {
                    if ((i == readOverlapLength - 1) && (! seen_snp)) {
                        non_vars_on_hmers++;
                    }
                }
            }

            String test = "Hello World";
        }
    }


    @Override
    public void closeTool() {
        // Print extra info
        outputStream.println("Hmers without SNPs: " + "\t" + non_vars_on_hmers);
        outputStream.println("SNPs on Edges: " + "\t" + snps_on_hmers_edge);

        // Print header
        outputStream.println("Hmer Position" + "\t" + "AtoC_PB" + "\t" + "AtoC_NP" + "\t" + "AtoG_PB" + "\t" + "AtoG_NP" + "\t" + "AtoT_PB" + "\t" + "AtoT_NP" + "\t" +
                "CtoA_PB" + "\t" + "CtoA_NP" + "\t" + "CtoG_PB" + "\t" + "CtoG_NP" + "\t" + "CtoT_PB" + "\t" + "CtoT_NP" + "\t" +
                "GtoA_PB" + "\t" + "GtoA_NP" + "\t" + "GtoC_PB" + "\t" + "GtoC_NP" + "\t" + "GtoT_PB" + "\t" + "GtoT_NP" + "\t" +
                "TtoA_PB" + "\t" + "TtoA_NP" + "\t" + "TtoC_PB" + "\t" + "TtoC_NP" + "\t" + "TtoG_PB" + "\t" + "TtoG_NP");

        // Print data to file
        for (int i = 0; i < 32; i++) {
            outputStream.print(i + "\t");
            outputStream.print(AtoC[i][0] + "\t" + AtoC[i][1] + "\t" + AtoG[i][0] + "\t" + AtoG[i][1] + "\t" + AtoT[i][0] + "\t" + AtoT[i][1] + "\t" +
                    CtoA[i][0] + "\t" + CtoA[i][1] + "\t" + CtoG[i][0] + "\t" + CtoG[i][1] + "\t" + CtoT[i][0] + "\t" + CtoT[i][1] + "\t" +
                    GtoA[i][0] + "\t" + GtoA[i][1] + "\t" + GtoC[i][0] + "\t" + GtoC[i][1] + "\t" + GtoT[i][0] + "\t" + GtoT[i][1] + "\t" +
                    TtoA[i][0] + "\t" + TtoA[i][1] + "\t" + TtoC[i][0] + "\t" + TtoC[i][1] + "\t" + TtoG[i][0] + "\t" + TtoG[i][1] + "\n");
        }

        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
