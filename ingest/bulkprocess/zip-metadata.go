package bulkprocess

import (
	"fmt"
	"path/filepath"
	"strings"
)

type ZipMetadata struct {
	SampleID string
	FieldID  string
	Instance string
	Index    string
}

func zipPathToMetadata(path string) (ZipMetadata, error) {
	filename := filepath.Base(path)

	// Remove .zip
	name := strings.Split(filename, ".")[0]

	data := strings.Split(name, "_")

	if len(data) != 4 {
		return ZipMetadata{}, fmt.Errorf("Expected filename to be of format sampleID_fieldID_instance_index.zip, but found %d parts instead of 4", len(data))
	}

	return ZipMetadata{data[0], data[1], data[2], data[3]}, nil
}
