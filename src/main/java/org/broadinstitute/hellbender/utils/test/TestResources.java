package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Class containing helper methods and resolved directories for test resources.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class TestResources {

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";

    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";
    public static final String publicTestDirRoot = publicTestDir.replace(publicTestDirRelative, "");

    public static final String packageRootTestDir = publicTestDir + "org/broadinstitute/hellbender/";
    public static final String toolsTestDir = packageRootTestDir + "tools/";

    public static final String GCS_GATK_TEST_RESOURCES = "gs://hellbender/test/resources/";

    public static final String GCS_b37_REFERENCE_2BIT = GCS_GATK_TEST_RESOURCES + "benchmark/human_g1k_v37.2bit";
    public static final String GCS_b37_CHR20_21_REFERENCE_2BIT = GCS_GATK_TEST_RESOURCES + "human_g1k_v37.20.21.2bit";

    /**
     * LARGE FILES FOR TESTING (MANAGED BY GIT LFS)
     */
    public static final String largeFileTestDir = new File(publicTestDir, "large").getAbsolutePath() + "/";

    // All of chromosomes 20 and 21 from the b37 reference
    public static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";

    public static final String b37_2bit_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.2bit";

    // All of chromosomes 20 and 21 from the b38 reference
    public static final String b38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";

    // ~600,000 reads from chromosomes 20 and 21 of an NA12878 WGS bam aligned to b37, plus ~50,000 unmapped reads
    public static final String NA12878_20_21_WGS_bam = largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";

    // Variants from a DBSNP 138 VCF overlapping the reads in NA12878_20_21_WGS_bam
    public static final String dbsnp_138_b37_20_21_vcf = largeFileTestDir + "dbsnp_138.b37.20.21.vcf";

    // Variants from a DBSNP 138 VCF form the first 65Mb of chr1
    public static final String dbsnp_138_b37_1_65M_vcf = largeFileTestDir + "dbsnp_138.b37.1.1-65M.vcf";

    public static final String WGS_B37_CH20_1M_1M1K_BAM = "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
    public static final String DBSNP_138_B37_CH20_1M_1M1K_VCF = "dbsnp_138.b37.excluding_sites_after_129.ch20.1m-1m1k.vcf";

    /**
     * END OF LARGE FILES FOR TESTING
     */

    public static final String NA12878_chr17_1k_BAM = publicTestDir + "NA12878.chr17_69k_70k.dictFix.bam";
    public static final String NA12878_chr17_1k_CRAM = publicTestDir + "NA12878.chr17_69k_70k.dictFix.cram";
    public static final String v37_chr17_1Mb_Reference = publicTestDir + "human_g1k_v37.chr17_1Mb.fasta";

    public static final String hg19_chr1_1M_Reference = publicTestDir + "Homo_sapiens_assembly19_chr1_1M.fasta";
    public static final String hg19_chr1_1M_dbSNP = publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";

    // the following file has been modified such that the first chromosome length is 1M; this is sometimes
    // required due to sequence dictionary validation, since a reference FASTA with only 1M bases is used
    public static final String hg19_chr1_1M_dbSNP_modified = publicTestDir + "HSA19.dbsnp135.chr1_1M.exome_intervals.modified.vcf";

    public static final String hg19_chr1_1M_exampleVCF = publicTestDir + "joint_calling.chr1_1M.1kg_samples.10samples.noINFO.vcf";
    public static final String hg19MiniReference = publicTestDir + "hg19mini.fasta";
    // Micro reference is the same as hg19mini, but contains only chromosomes 1 and 2
    public static final String hg19MicroReference = publicTestDir + "hg19micro.fasta";

    public static final String exampleFASTA = publicTestDir + "exampleFASTA.fasta";
    public static final String exampleReference = hg19MiniReference;
    public static final String hg19MiniIntervalFile = publicTestDir + "hg19mini.interval_list";


    private static CachingIndexedFastaSequenceFile hg19ReferenceReader;
    private static GenomeLocParser hg19GenomeLocParser;
    // used to seed the genome loc parser with a sequence dictionary
    public static SAMFileHeader hg19Header;


    /** Gets the genome location parser, initialized only once. */
    public static synchronized GenomeLocParser getGenomeLocParser() {
        if (hg19GenomeLocParser == null) {
            try {
                hg19ReferenceReader = new CachingIndexedFastaSequenceFile(new File(hg19MiniReference));
                hg19Header = new SAMFileHeader();
                hg19Header.setSequenceDictionary(hg19ReferenceReader.getSequenceDictionary());
                hg19GenomeLocParser = new GenomeLocParser(hg19ReferenceReader);
            } catch (FileNotFoundException e) {
                Assert.fail("Required file for test not found: " + hg19MiniReference);
            }
        }
        return hg19GenomeLocParser;
    }

    /** Used to seed the genome loc parser with a sequence dictionary. */
    public static synchronized SAMFileHeader getHg19Header() {
        getGenomeLocParser();
        return hg19Header;
    }

    /** Gets the HG19 reference reader, initialized only once. */
    public static synchronized CachingIndexedFastaSequenceFile getHg19ReferenceReader() {
        getGenomeLocParser();
        return hg19ReferenceReader;
    }

    public static List<GenomeLoc> intervalStringsToGenomeLocs(String... intervals) {
        return intervalStringsToGenomeLocs(Arrays.asList(intervals));
    }

    public static List<GenomeLoc> intervalStringsToGenomeLocs( List<String> intervals ) {
        List<GenomeLoc> locs = new ArrayList<>();
        for (String interval: intervals)
            locs.add(getGenomeLocParser().parseGenomeLoc(interval));
        return Collections.unmodifiableList(locs);
    }
}
