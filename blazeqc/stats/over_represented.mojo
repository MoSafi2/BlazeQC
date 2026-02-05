"""Over-represented sequence (split from stats_.mojo)."""

struct OverRepresentedSequence(Copyable, Movable):
    var seq: String
    var count: Int
    var percentage: Float64
    var hit: String

    fn __init__(
        out self,
        seq: String,
        count: Int,
        percentage: Float64,
        hit: String,
    ):
        self.seq = seq
        self.count = count
        self.percentage = percentage
        self.hit = hit

    fn __init__(
        out self,
        seq: StringLiteral,
        count: Int,
        percentage: Float64,
        hit: StringLiteral,
    ):
        self.seq = String(seq)
        self.count = count
        self.percentage = percentage
        self.hit = String(hit)
