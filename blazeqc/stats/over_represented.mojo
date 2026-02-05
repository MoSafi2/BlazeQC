"""Over-represented sequence (split from stats_.mojo)."""

@value
struct OverRepresentedSequence:
    var seq: String
    var count: Int
    var percentage: Float64
    var hit: String
