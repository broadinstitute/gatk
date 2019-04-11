package main

import (
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"path/filepath"
	"strings"
)

func main() {
	// Need:
	// - Folder containing filenames with filename containing the lookup, and text of LVEF
	// - LVEF file suffix.
	// - Filename with sample ID

	var LVEFFolder, SuffixValue, ManifestPath string
	var phenotype string
	flag.StringVar(&LVEFFolder, "lvef_folder", "", "Folder containing files that have phenotype values")
	flag.StringVar(&SuffixValue, "suffix_value", ".overlay", "Will split and keep only filename parts before this when extracting the dicom name")
	flag.StringVar(&ManifestPath, "manifest_path", "", "Path to the manifest file that will permit lookups from dicom name to sample_id")
	flag.StringVar(&phenotype, "phenotype", "", "Name of the column to assign for your phenotype (e.g., lvef)")

	flag.Parse()

	if LVEFFolder == "" || ManifestPath == "" || phenotype == "" {
		flag.PrintDefaults()
		return
	}

	lvefMap, err := BuildLVEFMap(LVEFFolder, SuffixValue)
	if err != nil {
		log.Fatalln(err)
	}

	manifestMap, err := BuildManifestMap(ManifestPath)
	if err != nil {
		log.Fatalln(err)
	}

	fmt.Printf("sample_id\t%s\n", phenotype)
	for dicomname, lvef := range lvefMap {
		if sampleID, exists := manifestMap[dicomname]; exists {
			fmt.Printf("%s\t%s\n", sampleID, lvef)
		}
	}
}

func BuildManifestMap(ManifestPath string) (map[string]string, error) {
	// dicom file name => sample_id
	ManifestMap := make(map[string]string)

	dicomCol := -1
	sampleIDCol := -1

	fileBytes, err := ioutil.ReadFile(ManifestPath)
	if err != nil {
		return nil, err
	}

	byteReader := bytes.NewReader(fileBytes)

	rdr := csv.NewReader(byteReader)
	rdr.Comma = '\t'
	for i := 0; ; i++ {
		line, err := rdr.Read()
		if err == io.EOF {
			break
		} else if err != nil {
			return nil, err
		}

		if i == 0 {
			for colID, col := range line {
				if col == "dicom_file" {
					dicomCol = colID
				} else if col == "sample_id" {
					sampleIDCol = colID
				}
			}

			if dicomCol == -1 || sampleIDCol == -1 {
				return nil, fmt.Errorf("Could not locate dicom name or sample_id from header. Headers are required for the manifest.")
			}
		}

		ManifestMap[line[dicomCol]] = line[sampleIDCol]
	}

	return ManifestMap, nil
}

func BuildLVEFMap(LVEFFolder, SuffixValue string) (map[string]string, error) {
	// dicom file name => string
	LVEFMap := make(map[string]string)

	files, err := ioutil.ReadDir(LVEFFolder)
	if err != nil {
		return nil, err
	}

	// Read the LVEFs, build the map of filename => LVEF
	for _, fileInfo := range files {
		if !strings.HasSuffix(fileInfo.Name(), ".txt") || !strings.Contains(fileInfo.Name(), SuffixValue) {
			continue
		}

		pieces := strings.Split(fileInfo.Name(), SuffixValue)
		if len(pieces) < 1 {
			continue
		}

		fileBytes, err := ioutil.ReadFile(filepath.Join(LVEFFolder, fileInfo.Name()))
		if err != nil {
			return nil, err
		}
		lvef := string(fileBytes)

		LVEFMap[pieces[0]] = strings.TrimSpace(lvef)
	}

	return LVEFMap, nil
}
