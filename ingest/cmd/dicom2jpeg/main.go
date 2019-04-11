package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"strconv"

	"image"
	"image/color"
	"image/jpeg"
	"math"
	"os"
	"sync"

	"github.com/gradienthealth/dicom"
	"github.com/gradienthealth/dicom/dicomlog"
	"github.com/gradienthealth/dicom/dicomtag"
)

// See 'batcher' for a process that can convert dicoms to jpegs in bulk fashion

var (
	printMetadata = flag.Bool("print-metadata", false, "Print image metadata")
	extractImages = flag.Bool("extract-images", false, "Extract images into separate files")
	verbose       = flag.Bool("verbose", false, "Activate high verbosity log operation")
)

// FrameBufferSize represents the size of the *Frame buffered channel for streaming calls
const FrameBufferSize = 100

func main() {
	// Update usage docs
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n%s <dicom file> [flags]\n", os.Args[0], os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()
	if len(flag.Args()) == 0 {
		flag.Usage()
		os.Exit(1)
	}
	if *verbose {
		dicomlog.SetLevel(math.MaxInt32)
	}
	path := flag.Arg(0)

	var parsedData *dicom.DataSet

	// Non-streaming parsing:
	p, err := dicom.NewParserFromFile(path, nil)
	if err != nil {
		log.Panic("error creating new parser", err)
	}
	parsedData, err = p.Parse(dicom.ParseOptions{
		DropPixelData: !*extractImages,
		ReturnTags: []dicomtag.Tag{
			// The overlay data
			{Group: 0x6000, Element: 0x3000},

			// Number of overlay rows
			{Group: 0x6000, Element: 0x0010},

			// Number of overlay columns
			{Group: 0x6000, Element: 0x0011},

			// Number of overlay frames
			{Group: 0x6000, Element: 0x0015},
			{Group: 0x6000, Element: 0x0014},

			// Bits allocated in the overlay
			{Group: 0x6000, Element: 0x0100},

			// Other metadata

			// Patient orientation
			{Group: 0x0020, Element: 0x0020},

			// Image Position
			{Group: 0x0020, Element: 0x0032},

			// Image Orientation
			{Group: 0x0020, Element: 0x0037},

			// Slice location
			{Group: 0x0020, Element: 0x0041},

			// Instance number
			{Group: 0x0020, Element: 0x0013},

			// Patient position
			{Group: 0x0018, Element: 0x5100},
		},
	})
	if parsedData == nil || err != nil {
		log.Panicf("Error reading %s: %v", path, err)
	}
	if *extractImages {
		for _, elem := range parsedData.Elements {
			if elem.Tag == dicomtag.PixelData {
				data := elem.Value[0].(dicom.PixelDataInfo)

				var wg sync.WaitGroup
				for frameIndex, frame := range data.Frames {
					wg.Add(1)
					go generateImage(&frame, frameIndex, &wg)
				}
				wg.Wait()

			}
		}
	}

	// Print Metadata from parsedData if needed
	if *printMetadata {
		log.Println(parsedData)
		var width, height uint16
		var instanceNumber string
		var patientX, patientY, patientZ float64
		hasOverlay := false
		for _, elem := range parsedData.Elements {
			if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x3000}) == 0 {
				hasOverlay = true
			}

			if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0010}) == 0 {
				width = elem.Value[0].(uint16)
			}

			if elem.Tag.Compare(dicomtag.Tag{Group: 0x6000, Element: 0x0011}) == 0 {
				height = elem.Value[0].(uint16)
			}

			if elem.Tag.Compare(dicomtag.Tag{Group: 0x0020, Element: 0x0013}) == 0 {
				instanceNumber = elem.Value[0].(string)
			}

			if elem.Tag.Compare(dicomtag.Tag{Group: 0x0020, Element: 0x0032}) == 0 {
				patientX, err = strconv.ParseFloat(elem.Value[0].(string), 32)
				if err != nil {
					continue
				}
				patientY, err = strconv.ParseFloat(elem.Value[1].(string), 32)
				if err != nil {
					continue
				}
				patientZ, err = strconv.ParseFloat(elem.Value[2].(string), 32)
				if err != nil {
					continue
				}
			}

			fmt.Printf("%s: VR %s: %v\n", elem.Tag, elem.VR, elem.String())
		}

		log.Printf("This is image %s. Overlay? %t. Overlay width %d x height %d.", instanceNumber, hasOverlay, width, height)
		log.Printf("Patient orientation:%.2f %.2f %.2f\n", patientX, patientY, patientZ)

	}

	log.Println("Complete.")
}

func generateImage(frame *dicom.Frame, frameIndex int, wg *sync.WaitGroup) {
	if frame.IsEncapsulated {
		go generateEncapsulatedImage(frame.EncapsulatedData, frameIndex, wg)
	} else {
		go generateNativeImage(frame.NativeData, frameIndex, wg)
	}
}

func generateEncapsulatedImage(frame dicom.EncapsulatedFrame, frameIndex int, wg *sync.WaitGroup) {
	defer wg.Done()
	path := fmt.Sprintf("image_%d.jpg", frameIndex) // TODO: figure out the image format
	ioutil.WriteFile(path, frame.Data, 0644)
	log.Printf("%s: %d bytes\n", path, len(frame.Data))
}

func generateNativeImage(frame dicom.NativeFrame, frameIndex int, wg *sync.WaitGroup) {
	defer wg.Done()
	i := image.NewGray16(image.Rect(0, 0, frame.Cols, frame.Rows))
	for j := 0; j < len(frame.Data); j++ {
		i.SetGray16(j%frame.Cols, j/frame.Rows, color.Gray16{Y: uint16(frame.Data[j][0])}) // for now, assume we're not overflowing uint16, assume gray image
	}
	name := fmt.Sprintf("image_%d.jpg", frameIndex)
	f, err := os.Create(name)
	if err != nil {
		fmt.Printf("Error while creating file: %s", err.Error())
	}
	jpeg.Encode(f, i, &jpeg.Options{Quality: 100})
	log.Printf("%s written \n", name)
}
