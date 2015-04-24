package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

//Illumina Dir Test Data
public final class BinTdUtil {
    public static final String ltStr(final int lane, final int tile) {
        return "s_" + lane + "_" + tile;
    }

    public static final byte A = (byte) 65;
    public static final byte C = (byte) 67;
    public static final byte G = (byte) 71;
    public static final byte T = (byte) 84;
    public static final byte P = (byte) 46; //dot
    public static final Map<String, List<ClusterData>> goldData = new HashMap<String, List<ClusterData>>();
    public static final Map<String, List<Integer>> goldIndices = new HashMap<String, List<Integer>>();
    public static final Map<String, Integer> goldSizes = new HashMap<String, Integer>();

    static {
        int lane = 1;
        int tile = 1101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 2, 10, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
                        makeCd(lane, tile, 1140, 2120, true,
                                new byte[]{P, C, C, C, C, A, A, C, A, T, T, C, T, A, A, T, T, A, T, G, C, C, T, C, A, C, A, A, C, T, C, T, C, T, T, T, T, T, T, T, T, T, T, T, T, T, A, A, C, T, T, T, G, C, A, A, A, T},
                                new byte[]{2, 16, 25, 33, 35, 37, 37, 35, 39, 37, 37, 35, 37, 40, 41, 41, 41, 40, 40, 41, 40, 40, 40, 40, 40, 31, 31, 31, 35, 35, 37, 35, 37, 31, 31, 31, 35, 35, 35, 35, 35, 39, 39, 39, 39, 37, 33, 31, 24, 37, 39, 40, 31, 33, 37, 39, 31, 31},
                                "CAACTCTC"),
                        makeCd(lane, tile, 1047, 2122, false,
                                new byte[]{P,C,T,A,A,P,G,P,A,C,T,P,T,G,P,G,T,G,T,G,C,P,P,P,P,P,P,P,A,P,P,P,P,P,P,T,C,A,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,15,26,31,31,2,19,2,18,31,31,2,18,31,2,17,27,31,31,31,31,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1069, 2159, true,
                                new byte[]{T,C,C,C,T,T,A,C,C,A,T,C,A,A,A,T,C,A,A,T,T,G,P,C,C,G,T,C,C,A,C,A,G,G,A,C,G,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{34,34,34,37,37,37,37,37,39,39,39,39,39,41,41,41,41,41,41,41,41,41,2,18,32,31,33,33,37,37,37,37,37,27,27,27,31,30,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                "GTCCACAG"),
                        makeCd(lane, tile, 1175, 2197, true,
                                new byte[]{C,C,C,C,T,G,A,G,G,A,C,A,C,C,A,T,C,C,C,A,C,T,C,C,A,C,C,A,A,C,A,T,T,A,A,G,A,G,C,T,G,G,G,G,A,A,C,A,T,C,C,A,G,A,A,A,G,G},
                                new byte[]{34,34,34,37,37,37,37,37,39,39,39,39,39,41,41,41,41,41,41,41,41,41,41,41,41,34,34,34,37,37,37,37,37,33,34,31,37,37,37,37,37,39,39,39,39,39,41,41,41,41,41,41,41,41,41,41,41,41},
                                "CCAACATT"),
                        makeCd(lane, tile, 1048, 2197, false,
                                new byte[]{P,C,T,C,C,P,G,P,T,C,A,P,C,A,P,G,T,G,G,A,G,P,P,P,P,P,P,P,C,P,P,P,P,P,P,G,T,G,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,15,26,30,31,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null)
                )
        );
        goldSizes.put(ltStr(lane, tile), 60);

        tile = 1201;
        goldIndices.put(ltStr(lane, tile), makeList(0, 1, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
                makeCd(lane, tile, 1187, 2100, true,
                        new byte[]{P,G,C,G,G,T,A,A,T,T,C,C,A,G,C,T,C,C,A,A,T,A,G,C,G,T,A,T,C,T,G,C,C,A,A,A,A,A,A,G,A,G,C,C,C,G,C,A,T,T,G,C,C,G,A,G,A,C},
                        new byte[]{2,16,25,33,33,17,31,35,39,39,37,39,39,40,40,40,40,39,39,40,40,38,39,38,38,34,34,34,37,37,37,37,37,28,27,28,26,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                        "TATCTGCC"),
                makeCd(lane, tile, 1045, 2105, false,
                        new byte[]{P,T,A,A,A,G,A,G,A,A,A,T,C,A,A,G,A,A,T,A,C,T,A,T,T,C,T,G,T,A,A,T,C,P,T,T,T,T,T,T,T,T,T,T,P,P,T,T,T,T,T,T,T,T,T,T,T,T},
                        new byte[]{2,12,19,31,30,7,31,8,31,31,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,16,16,16,17,31,32,2,2,2,15,26,31,31,31,31,31,31,31,30,2,2,15,25,30,30,30,30,30,30,30,30,28,27},
                        "CTGTAATC"),
                makeCd(lane, tile, 1159, 2179, false,
                        new byte[]{G,T,T,A,G,C,A,C,A,G,A,T,A,T,T,G,G,A,T,G,A,G,T,G,A,A,A,A,A,A,A,A,A,T,T,T,T,T,T,T,T,T,A,T,T,T,T,T,C,T,A,A,A,T,A,C,T,T},
                        new byte[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,28,28,28,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                        null),
                makeCd(lane, tile, 1103, 2184, true,
                        new byte[]{G,T,A,A,G,A,A,C,T,A,C,C,C,T,G,G,G,T,C,C,C,C,G,T,G,T,T,G,T,C,T,A,T,A,G,A,A,G,T,T,T,C,A,G,A,A,T,T,G,T,G,G,C,C,C,C,A,T},
                        new byte[]{31,31,33,37,37,37,37,37,39,39,39,39,39,41,41,41,41,38,40,41,41,41,41,39,40,31,34,34,37,37,37,37,37,33,31,33,37,37,35,36,37,39,39,39,39,39,41,41,41,38,39,40,41,41,41,41,41,40},
                        "TTGTCTAT")
        ));
        goldSizes.put(ltStr(lane, tile), 60);

        tile = 2101;
        goldIndices.put(ltStr(lane, tile), makeList(7, 15, 16, 19));
        goldData.put(ltStr(lane, tile), makeList(
                makeCd(lane, tile, 1123, 2095, true,
                        new byte[]{P,T,G,G,A,C,A,A,C,A,T,G,T,T,C,G,A,G,A,G,C,T,A,C,A,C,A,G,C,G,G,T,A,T,C,C,G,C,C,T,C,C,A,G,C,T,T,C,A,G,C,T,T,C,T,C,C,T},
                        new byte[]{2,16,28,33,33,35,35,35,37,37,37,37,35,38,37,38,40,38,30,37,26,39,39,37,40,31,30,31,35,35,37,31,31,31,31,31,37,35,35,37,37,39,39,39,39,39,41,39,38,38,41,40,41,41,41,36,39,39},
                        "CAGCGGTA"),
               makeCd(lane, tile, 1162, 2139, true,
                       new byte[]{A,G,A,G,G,T,G,A,A,A,T,T,C,T,T,G,G,A,C,C,G,G,C,G,C,T,G,C,T,G,C,T,G,A,T,C,G,T,T,T,A,T,G,G,T,C,G,G,A,A,C,T,A,C,G,A,C,G},
                       new byte[]{31,31,31,35,35,35,35,35,39,37,39,39,39,35,33,25,36,37,39,39,34,32,38,30,35,34,34,34,37,37,37,37,37,33,34,34,37,37,37,37,37,39,39,39,39,39,40,41,41,41,41,41,41,41,40,41,41,40},
                       "TGCTGCTG"),
                makeCd(lane, tile, 1013, 2146, true,
                        new byte[]{P,A,C,A,C,T,G,C,T,G,C,A,G,A,T,G,A,C,A,A,G,C,A,G,C,C,T,A,T,G,C,G,T,P,P,P,P,C,G,C,T,A,G,A,A,C,C,A,A,C,T,T,A,T,T,C,A,T},
                        new byte[]{2,19,33,35,37,37,37,37,39,39,39,39,39,41,41,41,41,41,41,41,41,41,41,41,41,34,34,34,37,37,37,37,37,2,2,2,2,17,19,28,30,31,31,30,31,30,31,31,30,31,31,31,31,31,31,30,31,31},
                        "CTATGCGT"),
                makeCd(lane, tile, 1245, 2154, true,
                        new byte[]{T,C,G,T,T,A,A,G,T,A,T,A,T,T,C,T,T,A,G,G,T,A,T,T,T,C,T,G,T,A,A,T,C,A,C,C,A,A,T,C,A,G,T,A,G,C,A,C,C,A,C,T,A,T,A,C,A,C},
                        new byte[]{34,34,34,37,37,35,37,37,37,39,37,39,39,40,40,41,41,41,41,41,37,41,41,41,40,31,34,34,37,37,37,37,37,34,34,34,37,37,37,37,37,39,39,39,39,39,41,41,41,41,41,41,40,41,41,41,41,41},
                        "CTGTAATC")
        ));
        goldSizes.put(ltStr(lane, tile), 60);

        tile = 11101;
        goldIndices.put(ltStr(lane, tile), makeList(0, 2, 10, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
                        makeCd(lane, tile, 1140, 2120, true,
                                new byte[]{P,A,A,C,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1047, 2122, false,
                                new byte[]{P,A,A,G,A,C,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,27,37,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1069, 2159, true,
                                new byte[]{P,C,T,T,G,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,37,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1175, 2197, true,
                                new byte[]{P,A,A,A,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1048, 2197, false,
                                new byte[]{P,A,A,G,A,C,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,27,22,32,32,37,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null)
                )
        );
        goldSizes.put(ltStr(lane, tile), 341292);

        tile = 11102;
        goldIndices.put(ltStr(lane, tile), makeList(0, 2, 10, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
                        makeCd(lane, tile, 1140, 2120, true,
                                new byte[]{P,A,A,G,A,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1047, 2122, false,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1069, 2159, true,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1175, 2197, true,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1048, 2197, false,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null)
                )
        );
        goldSizes.put(ltStr(lane, tile), 366884);

        tile = 11103;
        goldIndices.put(ltStr(lane, tile), makeList(0, 2, 10, 18, 19));
        goldData.put(ltStr(lane, tile), makeList(
                        makeCd(lane, tile, 1140, 2120, true,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1047, 2122, false,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1069, 2159, true,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1175, 2197, true,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null),
                        makeCd(lane, tile, 1048, 2197, false,
                                new byte[]{P,G,C,T,T,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P,P},
                                new byte[]{2,32,32,32,32,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
                                null)
                )
        );
        goldSizes.put(ltStr(lane, tile), 336434);
    }

    public static Map<Integer, ClusterData> clusterData(final int lane, final List<Integer> tiles, final String readStructure, final IlluminaDataType... dataTypes) {
        final List<Integer> sortedTiles = new ArrayList<Integer>(tiles);
        Collections.sort(sortedTiles);

        final Map<Integer, ClusterData> data = new HashMap<Integer, ClusterData>();
        int offset = 0;
        for (final int tile : sortedTiles) {
            final String key = ltStr(lane, tile);
            final List<ClusterData> cds = goldData.get(key);
            final List<Integer> readNos = goldIndices.get(key);
            final int size = goldSizes.get(key);

            for (int i = 0; i < cds.size(); i++) {
                data.put(offset + readNos.get(i), selectiveCopyCd(cds.get(i), readStructure, dataTypes));
            }

            offset += size;
        }
        return data;
    }

    public static ReadData[] copyReadData(final ReadStructure rs, final IlluminaDataType[] dts, final ClusterData toCopy) {
        boolean doBases = false;
        boolean doQuals = false;

        for (final IlluminaDataType dt : dts) {
            switch (dt) {
                case BaseCalls:
                    doBases = true;
                    break;

                case QualityScores:
                    doQuals = true;
                    break;
            }
        }

        if (!doBases && !doQuals)
            return null;

        final ReadData rdToCopy = toCopy.getRead(0); //Only gonna be one read in this
        final ReadData[] rds = new ReadData[rs.nonSkips.length()];

        int index = 0;
        int baseIndex = 0;
        for (int i = 0; i < rs.descriptors.size(); i++) {
            final ReadDescriptor readDesc = rs.descriptors.get(i);

            if (readDesc.type != ReadType.S) {
                final ReadData curRead = new ReadData(readDesc.type);
                if (doBases) {
                    final byte[] bases = Arrays.copyOfRange(rdToCopy.getBases(), baseIndex, baseIndex + readDesc.length);
                    curRead.setBases(bases);
                }

                if (doQuals) {
                    final byte[] quals = Arrays.copyOfRange(rdToCopy.getQualities(), baseIndex, baseIndex + readDesc.length);
                    curRead.setQualities(quals);
                }

                rds[index++] = curRead;
            }

            baseIndex += readDesc.length;
        }

        return rds;
    }

    private static FourChannelIntensityData copyIntensities(final FourChannelIntensityData toCopy, final int start, final int length) {
        final FourChannelIntensityData fcid = new FourChannelIntensityData(length);

        System.arraycopy(toCopy.getA(), start, fcid.getA(), 0, length);
        System.arraycopy(toCopy.getC(), start, fcid.getC(), 0, length);
        System.arraycopy(toCopy.getG(), start, fcid.getG(), 0, length);
        System.arraycopy(toCopy.getT(), start, fcid.getT(), 0, length);
        return fcid;
    }

    public static ClusterData selectiveCopyCd(final ClusterData toCopy, final String readStructure, final IlluminaDataType... dataTypes) {
        final ReadStructure rs = new ReadStructure(readStructure);
        final ReadData[] rd = copyReadData(rs, dataTypes, toCopy);
        final ClusterData cd = new ClusterData(rd);
        cd.setTile(toCopy.getTile());
        cd.setLane(toCopy.getLane());

        for (final IlluminaDataType idt : dataTypes) {
            switch (idt) {
                case Position:
                    cd.setX(toCopy.getX());
                    cd.setY(toCopy.getY());
                    break;

                case PF:
                    cd.setPf(toCopy.isPf());
                    break;

                case Barcodes:
                    cd.setMatchedBarcode(toCopy.getMatchedBarcode());
                    break;

                default:
                    break;
            }
        }

        return cd;
    }

    public static ClusterData makeCd(final int lane, final int tile, final int xCoord, final int yCoord, final boolean pf, final byte[] bases, final byte[] qualities, final String matchedBarcode) {
        final ReadData rd = new ReadData();
        rd.setBases(Arrays.copyOf(bases, bases.length));
        rd.setQualities(Arrays.copyOf(qualities, bases.length));
        rd.setReadType(ReadType.T); //This will be ignored, as the cluster will be chopped up by ReadStructure

        final ClusterData cd = new ClusterData(rd);
        cd.setLane(lane);
        cd.setTile(tile);
        cd.setX(xCoord);
        cd.setY(yCoord);
        cd.setPf(pf);
        cd.setMatchedBarcode(matchedBarcode);

        return cd;
    }

}
