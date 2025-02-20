package org.broadinstitute.hellbender.tools.sv;

import java.io.*;
import java.util.*;

public class Pedigree {

    public enum Sex {
        MALE("1"), FEMALE("2"), UNKNOWN("other");

        /// stores PED format compatible format for to/from string conversions
        private final String pedEncode;
        Sex(String code) {
            this.pedEncode = code;
        }

        public static Sex fromString(String sex) {
            return switch (sex) {
                case "1" -> MALE;
                case "2" -> FEMALE;
                default -> UNKNOWN;
            };
        }

        @Override
        public String toString() {
            return pedEncode;
        }
    }

    public enum Phenotype {
        MISSING("0"), UNAFFECTED("1"), AFFECTED("2");

        /// stores PED format compatible format for to/from string conversions
        private final String pedEncode;
        Phenotype(String code) {
            this.pedEncode = code;
        }

        public static Phenotype fromString(String affected) {
            return switch (affected) {
                case "-9", "0" -> MISSING;
                case "1" -> UNAFFECTED;
                case "2" -> AFFECTED;
                default -> throw new IllegalArgumentException("Invalid affected status: " + affected);
            };
        }

        @Override
        public String toString() {
            return pedEncode;
        }
    }

    public record Relationship(
            String familyId,
            String individualId,
            String fatherId,
            String motherId,
            Sex sex,
            Phenotype phenotype
    ) {
        @Override
        public String toString() {
            return toString("\t");
        }

        public String toString(String delimiter) {
            return String.join(delimiter,
                    familyId,
                    individualId,
                    fatherId,
                    motherId,
                    sex.toString(),
                    phenotype.toString()
            );
        }
    }

    private final Map<String, Relationship> trios = new HashMap<>();

    public Pedigree() { }

    public Pedigree(String pedFilename) throws Exception
    {
        var pedFile = new File(pedFilename);
        if(!pedFile.exists()) {
            throw new FileNotFoundException(pedFilename + " not found");
        }
        else if(pedFile.isDirectory()) {
            throw new IllegalArgumentException("Invalid filename, " + pedFilename + " is a directory, not a file");
        }

        try (var reader = new BufferedReader(new FileReader(pedFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                var parts = line.split("\t");
                var individual = parts[1];
                trios.put(
                    individual,
                    new Relationship(
                        parts[0],
                        individual,
                        parts[2],
                        parts[3],
                        Sex.fromString(parts[4]),
                        Phenotype.fromString(parts[5])));
            }
        } catch (IOException e) {
            throw new IOException(pedFilename + " does not match the expected format");
        }
    }

    public void addIndividual(Relationship relationship) {
        trios.put(relationship.individualId, relationship);
    }

    public Integer getTrioCount() {
        return trios.size();
    }

    public Relationship getTrio(String individualId) {
        return trios.get(individualId);
    }

    public void toPed(String pedFilename) throws IOException {
        try (var writer = new BufferedWriter(new FileWriter(pedFilename))) {
            for (var relationship : trios.entrySet()) {
                writer.write(relationship.getValue().toString());
                writer.newLine();
            }
        } catch (IOException e) {
            throw new IOException("Error writing to file " + pedFilename + e.getMessage());
        }
    }
}
