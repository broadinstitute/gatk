package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"

	"github.com/carbocation/genomisc"
)

const (
	// .sample file field columns
	ID_1 = iota
	ID_2
	missing
	sex
)

func main() {
	var (
		samplePath string
	)

	flag.StringVar(&samplePath, "sample", "", "genotyping .sample file for the UKBB")
	flag.Parse()

	if samplePath == "" {
		flag.PrintDefaults()
		os.Exit(1)
	}

	samplePath = genomisc.ExpandHome(samplePath)
	log.Printf("Importing %s\n", samplePath)

	// .sample file

	f, err := os.Open(samplePath)
	if err != nil {
		log.Fatalln(err)
	}
	defer f.Close()

	delim := genomisc.DetermineDelimiter(f)

	f.Seek(0, 0)
	fileCSV := csv.NewReader(f)
	fileCSV.Comma = delim

	// .sample files have 2 header rows that we will discard
	fileCSV.Read()
	fileCSV.Read()

	i := 0
	fmt.Printf("sample_id\tfile_row\n")
	for ; ; i++ {
		row, err := fileCSV.Read()
		if err != nil && err == io.EOF {
			break
		} else if err != nil {
			log.Fatalln(err)
		}

		fmt.Printf("%s\t%d\n", row[ID_1], i)
	}

	log.Println("Extracted", i, "records from the .sample file")
}
