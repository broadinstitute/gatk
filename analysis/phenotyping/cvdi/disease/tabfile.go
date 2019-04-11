package main

import (
	"encoding/csv"
	"errors"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type TabEntry struct {
	FieldID int
	Values  []string
	Exclude bool
}

// The UK Biobank keys up the ICD and OPCS codes without any decimals, so e.g.,
// K41.2 becomes K412. You can look up the true value using the coding table,
// but for simplicity we just strip the dots.
func (t TabEntry) FormattedValues() []string {
	out := make([]string, 0, len(t.Values))

	// HESIN is special-cased to exclude "." (So K41.2 becomes K412)
	if IsHesin(t.FieldID) {
		for _, v := range t.Values {
			out = append(out, strings.Replace(v, ".", "", -1))
		}
	} else {
		out = append(out, t.Values...)
	}

	// Every field will get leading and trailing spaces trimmed
	for i, v := range out {
		out[i] = strings.TrimSpace(v)
	}

	return out
}

type TabFile struct {
	Include struct {
		Hesin    []TabEntry
		Special  []TabEntry
		Standard []TabEntry
	}

	Exclude struct {
		Hesin    []TabEntry
		Special  []TabEntry
		Standard []TabEntry
	}
}

func NewTabFile() *TabFile {
	t := &TabFile{}

	t.Include.Hesin = make([]TabEntry, 0)
	t.Include.Special = make([]TabEntry, 0)
	t.Include.Standard = make([]TabEntry, 0)

	t.Exclude.Hesin = make([]TabEntry, 0)
	t.Exclude.Special = make([]TabEntry, 0)
	t.Exclude.Standard = make([]TabEntry, 0)

	return t
}

func (t *TabFile) AllStandardFields() []int {
	all := make([]int, 0)

	seen := make(map[int]struct{})

	for _, v := range append(t.Include.Standard, t.Exclude.Standard...) {
		if _, exists := seen[v.FieldID]; exists {
			continue
		}

		all = append(all, v.FieldID)
	}

	return all
}

func (t *TabFile) AllIncluded() []TabEntry {
	all := append(t.Include.Hesin, t.Include.Special...)
	return append(all, t.Include.Standard...)
}

func (t *TabFile) AllExcluded() []TabEntry {
	all := append(t.Exclude.Hesin, t.Exclude.Special...)
	return append(all, t.Exclude.Standard...)
}

func (t *TabFile) CheckSensibility() ([]string, error) {
	batchnames := []string{"Inclusion", "Exclusion"}
	batches := [][]TabEntry{t.AllIncluded(), t.AllExcluded()}

	missingFields := make(map[int]struct{})
	errs := make([]string, 0)

	for batchID, batch := range batches {
		icd9 := make(map[int]struct{})
		icd10 := make(map[int]struct{})
		opcs := make(map[int]struct{})

		for _, v := range batch {
			if _, exists := ICD9[v.FieldID]; exists {
				icd9[v.FieldID] = struct{}{}
			}
			if _, exists := ICD10[v.FieldID]; exists {
				icd10[v.FieldID] = struct{}{}
			}
			if _, exists := OPCS[v.FieldID]; exists {
				opcs[v.FieldID] = struct{}{}
			}
		}

		if len(icd9) > 0 && len(icd9) < len(ICD9) {
			for known := range ICD9 {
				if _, exists := icd9[known]; !exists {
					missingFields[known] = struct{}{}
				}
			}
			errs = append(errs, fmt.Sprintf("In the %s subset, you included ICD9 fields %v. Consider including all IC9 fields: %v", batchnames[batchID], icd9, ICD9))
		}

		if len(icd10) > 0 && len(icd10) < len(ICD10) {
			for known := range ICD10 {
				if _, exists := icd10[known]; !exists {
					missingFields[known] = struct{}{}
				}
			}
			errs = append(errs, fmt.Sprintf("In the %s subset, you included ICD10 fields %v. Consider including all IC10 fields: %v", batchnames[batchID], icd10, ICD10))
		}

		if len(opcs) > 0 && len(opcs) < len(OPCS) {
			for known := range OPCS {
				if _, exists := opcs[known]; !exists {
					missingFields[known] = struct{}{}
				}
			}
			errs = append(errs, fmt.Sprintf("In the %s subset, you included OPCS fields %v. Consider including all OPCS fields: %v", batchnames[batchID], opcs, OPCS))
		}
	}

	if len(errs) > 0 {
		missingSlice := make([]string, 0, len(missingFields))
		for missing := range missingFields {
			missingSlice = append(missingSlice, strconv.Itoa(missing))
		}

		return missingSlice, errors.New(strings.Join(errs, " | "))
	}

	return nil, nil
}

func ParseTabFile(tabPath string) (*TabFile, error) {
	f, err := os.Open(tabPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	fileCSV := csv.NewReader(f)
	fileCSV.Comma = '\t'

	recs, err := fileCSV.ReadAll()
	if err != nil {
		log.Fatalln(err)
	}

	output := NewTabFile()

	for i, row := range recs {
		if l := len(row); l != 3 {
			return nil, fmt.Errorf("Tabfile %s row %d had %d columns, expected 3", tabPath, i, l)
		}

		if i == 0 {
			// header
			continue
		}
		entry := TabEntry{
			Values:  strings.Split(row[1], ","),
			Exclude: row[2] == "1",
		}
		entry.FieldID, err = strconv.Atoi(row[0])
		if err != nil {
			return nil, err
		}

		// Assign to the right field type
		switch entry.Exclude {
		case true:
			if IsHesin(entry.FieldID) {
				output.Exclude.Hesin = append(output.Exclude.Hesin, entry)
			} else if IsSpecial(entry.FieldID) {
				output.Exclude.Special = append(output.Exclude.Special, entry)
			} else {
				output.Exclude.Standard = append(output.Exclude.Standard, entry)
			}
		default:
			if IsHesin(entry.FieldID) {
				output.Include.Hesin = append(output.Include.Hesin, entry)
			} else if IsSpecial(entry.FieldID) {
				output.Include.Special = append(output.Include.Special, entry)
			} else {
				output.Include.Standard = append(output.Include.Standard, entry)
			}
		}
	}

	return output, nil
}
