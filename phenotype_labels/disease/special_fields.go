package main

var (
	// These identify which FieldIDs should be special-cased: either to the
	// HESIN table or the table with fields that have a specially-known date.
	// All other fields can be queried against the main phenotype table, with an
	// assumed date of whatever makes the most sense for your purposes
	// (generally, I use the enrollment date, but you could use birthdate, etc).

	MaterializedHesin = map[int]struct{}{
		41210: struct{}{},
		41202: struct{}{},
		41204: struct{}{},
		40001: struct{}{},
		40002: struct{}{},
		41200: struct{}{},
		41203: struct{}{},
		41205: struct{}{},
	}

	MaterializedSpecial = map[int]struct{}{
		42013: struct{}{},
		42011: struct{}{},
		42009: struct{}{},
		42007: struct{}{},
		42001: struct{}{},
		20004: struct{}{},
		20002: struct{}{},
		20001: struct{}{},
	}
)

var (
	// These can be helpful to make sure that the user is including all fields
	// that use the same family of codes

	ICD9 = map[int]struct{}{
		41203: struct{}{},
		41205: struct{}{},
	}

	ICD10 = map[int]struct{}{
		41202: struct{}{},
		41204: struct{}{},
		40001: struct{}{},
		40002: struct{}{},
	}

	OPCS = map[int]struct{}{
		41200: struct{}{},
		41210: struct{}{},
	}
)

func IsHesin(fieldID int) bool {
	_, exists := MaterializedHesin[fieldID]

	return exists
}

func IsSpecial(fieldID int) bool {
	_, exists := MaterializedSpecial[fieldID]

	return exists
}
