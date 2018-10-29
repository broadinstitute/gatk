"""
Constants that must remain in sync with the companion StreamingProcessController Java
code in GATK. See StreamingToolConstants.java.

"""

"""
Command acknowledgement messages used to signal positive acknowledgement ('ack',
negative acknowledgement ('nck'), and negative acknowledgement with an accompanying
message ('nkm').
"""
_ackString = "ack"
_nackString = "nck"
_nkmString = "nkm"


"""
The length of a message written with a negative ack (nkm) must be 4 bytes long when
serialized as a string, and cannot have a value > 9999.
"""
_nckMessageLengthSerializedSize = 4
_nckMaxMessageLength = 9999
