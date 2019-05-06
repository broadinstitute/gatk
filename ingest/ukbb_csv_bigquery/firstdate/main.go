package main

func main() {
	// http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41253

	// Run and merge from 3 separate queries:
	//
	// 1) for ICD/opcode fields, look into hesin
	// 2) for special fields with dates, look at their dates
	// 3) for all other fields, use the enrollment date based on their array_idx

	// assign ICDs to their transformed FieldIDs - different for main and secondary

	// Then, the output:
	// SampleID FieldID Value Date

	// Then left join on censor
	// GROUP BY sample_id

	// Downstream:
	// Pass 1: Fetch censor data for everyone
	// Pass 2:
}
