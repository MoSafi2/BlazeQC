"""Over-represented sequence (split from stats_.mojo)."""

struct OverRepresentedSequence(Copyable, Movable):
    var seq: String
    var count: Int
    var percentage: Float64
    var hit: String
