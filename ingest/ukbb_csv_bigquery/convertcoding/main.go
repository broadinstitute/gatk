package main

import (
	"bytes"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"strings"
)

const (
	ExpectedRows = 3
)

func main() {
	var (
		codingPath string
	)

	flag.StringVar(&codingPath, "coding", "https://biobank.ctsu.ox.ac.uk/~bbdatan/Codings_Showcase.csv", "URL to CSV file with the UKBB data encodings")
	flag.Parse()

	if codingPath == "" {
		flag.PrintDefaults()
		log.Fatalln()
	}

	if err := ImportCoding(codingPath); err != nil {
		log.Fatalln(err)
	}
}

func ImportCoding(url string) error {
	log.Printf("Importing from %s\n", url)

	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	reader := csv.NewReader(resp.Body)
	reader.Comma = ','
	reader.LazyQuotes = true

	header := make([]string, 0)
	j := 0
	for ; ; j++ {
		row, err := reader.Read()
		if err != nil && err == io.EOF {
			resp.Body.Close()
			break
		} else if err != nil {
			buf := bytes.NewBuffer(nil)
			io.Copy(buf, resp.Body)
			if strings.Contains(buf.String(), "internal error") {
				log.Println("Coding File is not permitted to be downloaded from the UKBB")
				continue
			}
		}

		// Handle the header
		if j == 0 {
			log.Printf("Header (%d elements): %+v\n", len(row), row)
			header = append(header, row...)
			for k, v := range header {
				if v == "Coding" {
					header[k] = "coding_file_id"
				} else if v == "Value" {
					header[k] = "coding"
				} else if v == "Meaning" {
					header[k] = "meaning"
				}
			}

			if nCols := len(header); nCols != ExpectedRows {
				return fmt.Errorf("Expected a CSV with %d columns; got one with %d", ExpectedRows, nCols)
			}

			fmt.Println(strings.Join(header, "\t"))

			continue
		}

		// Handle the entries
		if len(row) == ExpectedRows {
			fmt.Println(strings.Join(row, "\t"))
		}
	}

	log.Println("Created coding file with", j, "entries")

	return nil
}
