from blazeseq.quality_schama import (
    QualitySchema,
    sanger_schema,
    illumina_1_3_schema,
    solexa_schema,
    illumina_1_5_schema,
    illumina_1_8_schema,
    generic_schema
)
from utils.variant import Variant

alias schema = Variant[String, QualitySchema]
alias read_header = ord("@")
alias quality_header = ord("+")
alias new_line = ord("\n")
alias carriage_return = ord("\r")


struct FastqRecord(Copyable, Movable, Writable 
#Sized, Stringable, Writable
):
    """Struct that represent a single FastaQ record."""

    var SeqHeader: String
    var SeqStr: String
    var QuHeader: String
    var QuStr: String
    var quality_schema: QualitySchema

    fn __init__(
        out self,
        SeqHeader: String,
        SeqStr: String,
        QuHeader: String,
        QuStr: String,
        quality_schema: schema = "generic",
    ) raises:
        self.SeqHeader = SeqHeader
        self.QuHeader = QuHeader
        self.SeqStr = SeqStr
        self.QuStr = QuStr

        if quality_schema.isa[String]():
            self.quality_schema = self._parse_schema(quality_schema[String])
        else:
            self.quality_schema = quality_schema[QualitySchema]

    @always_inline
    fn get_seq(self) -> StringSlice[__origin_of(self.SeqStr)]:
        return self.SeqStr.as_string_slice()

    @always_inline
    fn get_quality_string(self) -> StringSlice[__origin_of(self.QuStr)]:
        return self.QuStr.as_string_slice()

    @always_inline
    fn get_qulity_scores(
        self, quality_format: String
    ) -> List[UInt8]:
        var schema = self._parse_schema((quality_format))
        output = List[UInt8](capacity=len(self.QuStr))
        bytes = self.QuStr.as_bytes()
        for i in range(len(self.QuStr)):
            output[i] = bytes[i] - schema.OFFSET
        return output

    @always_inline
    fn get_qulity_scores(self, schema: QualitySchema) -> List[UInt8]:
        output = List[UInt8](capacity=len(self.QuStr))
        bytes = self.QuStr.as_bytes()
        for i in range(len(self.QuStr)):
            output[i] = bytes[i] - schema.OFFSET
        return output

    @always_inline
    fn get_qulity_scores(self, offset: UInt8) -> List[UInt8]:
        output = List[UInt8](capacity=len(self.QuStr))
        bytes = self.QuStr.as_bytes()
        for i in range(len(self.QuStr)):
            output[i] = bytes[i] - offset
        return output

    @always_inline
    fn get_header_string(self) -> StringSlice[__origin_of(self.SeqHeader)]:
        return self.SeqHeader.as_string_slice()

    @always_inline
    fn validate_record(self) raises:
        if self.SeqHeader.as_bytes()[0] != read_header:
            raise Error("Sequence Header is corrupt")

        if self.QuHeader.as_bytes()[0] != quality_header:
            raise Error("Quality Header is corrupt")

        if len(self.SeqStr) != len(self.QuStr):
            raise Error("Corrput Lengths")

        if len(self.QuHeader) > 1:
            if len(self.QuHeader) != len(self.SeqHeader):
                raise Error("Quality Header is corrupt")
            
            if not self.QuHeader.as_string_slice()[1:] != self.SeqHeader.as_string_slice()[1:]:
                raise Error("Quality Header is not the same as the Sequecing Header")

    @always_inline
    fn validate_quality_schema(self) raises:
        for i in range(len(self.QuStr)):
            if (
                self.QuStr.as_bytes()[i] > self.quality_schema.UPPER
                or self.QuStr.as_bytes()[i] < self.quality_schema.LOWER
            ):
                raise Error(
                    "Corrput quality score according to proivded schema"
                )

    @always_inline
    fn total_length(self) -> Int:
        return (
            len(self.QuHeader)
            + len(self.QuStr)
            + len(self.SeqHeader)
            + len(self.SeqStr)
        )


    @staticmethod
    @always_inline
    fn _parse_schema(quality_format: String) -> QualitySchema:
        var schema: QualitySchema

        if quality_format == "sanger":
            schema = sanger_schema
        elif quality_format == "solexa":
            schema = solexa_schema
        elif quality_format == "illumina_1.3":
            schema = illumina_1_3_schema
        elif quality_format == "illumina_1.5":
            schema = illumina_1_5_schema
        elif quality_format == "illumina_1.8":
            schema = illumina_1_8_schema
        elif quality_format == "generic":
            schema = generic_schema
        else:
            print(
                "Uknown quality schema please choose one of 'sanger', 'solexa',"
                " 'illumina_1.3', 'illumina_1.5' 'illumina_1.8', or 'generic'"
            )
            return generic_schema
        return schema

    # BUG: returns Smaller strings that expected.
    @always_inline
    fn __str__(self) -> String:
        return String.write(self)

    
    fn write_to[w: Writer](self, mut writer: w):
        writer.write(self.SeqHeader, "\n", self.SeqStr, "\n", self.QuHeader, "\n", self.QuStr, "\n")


    @always_inline
    fn __len__(self) -> Int:
        return len(self.SeqStr)

    # @always_inline
    # fn hash[bits: Int = 3, length: Int = 64 // bits](self) -> UInt64:
    #     """Hashes the first xx bp (if possible) into one 64bit. Max length is 64/nBits per bp.
    #     """

    #     @parameter
    #     if length < 32:
    #         return self._hash_packed(self.SeqStr.unsafe_ptr(), length)
    #     return self._hash_additive(self.SeqStr.unsafe_ptr(), length)

    # # Can be Vectorized
    # @staticmethod
    # @always_inline
    # fn _hash_packed[
    #     bits: Int = 3
    # ](bytes: UnsafePointer[Byte], length: Int) -> UInt64:
    #     """
    #     Hash the DNA strand to into 64bits unsigned number using xbit encoding.
    #     If the length of the bytes strand is longer than 64//bits bps, the hash is truncated.
    #     ----

    #     parameters:
    #     - bits (Int): the number of least significant bits used to hash a base pair. increased bit width reduces the number of bp that can be hashed.

    #     args:
    #     - bytes (UnsafePointer[Byte]): pointer the the basepair buffer.
    #     - length (Int): the length of the buffer to be hashed.
    #     """
    #     alias rnge: Int = 64 // bits
    #     alias width = simdwidthof[Byte]()
    #     var hash: UInt64 = 0
    #     var mask = (0b1 << bits) - 1
    #     for i in range(min(rnge, length)):
    #         # Mask for for first <n> significant bits, vectorized operation.
    #         var base_val = bytes[i] & mask
    #         hash = (hash << bits) | Int(base_val[i])
    #     return hash

    # # Change to a better hashing Algorithm
    # @staticmethod
    # @always_inline
    # fn _hash_additive[
    #     bits: Int = 3
    # ](bytes: UnsafePointer[UInt8], length: Int) -> UInt64:
    #     """Hashes longer DNA sequences . It hashes 16bps spans of the sequences and using 2 or 3 bit encoding and adds them to the hash.
    #     """
    #     constrained[
    #         bits <= 3, "Additive hashing can only hash up to 3bit resolution"
    #     ]()
    #     var full_hash: UInt64 = 0
    #     var mask = (0b1 << bits) - 1
    #     var rounds = align_down(length, 16)
    #     var rem = length % 16

    #     for round in range(rounds):
    #         var interim_hash: UInt64 = 0

    #         @parameter
    #         for i in range(16):
    #             var base_val = bytes[i + 16 * round] & mask
    #             interim_hash = interim_hash << bits | Int(base_val)
    #         full_hash = full_hash + interim_hash

    #     if rem > 0:
    #         var interim_hash: UInt64 = 0
    #         for i in range(rem):
    #             var base_val = bytes[i + 16 * rounds] & mask
    #             interim_hash = interim_hash << bits | Int(base_val)
    #         full_hash = full_hash + interim_hash

    #     return full_hash

    # @always_inline
    # fn __hash__(self) -> UInt:
    #     return Int(self.hash())

    # @always_inline
    # fn __eq__(self, other: Self) -> Bool:
    #     return self.__hash__() == other.__hash__()

    # fn __ne__(self, other: Self) -> Bool:
    #     return self.__hash__() != other.__hash__()
