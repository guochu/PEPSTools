abstract type AbstractMessage end


"""
	struct Message{_I, _O}

	on a 2D lattice, the input message i is defined from left to right, or from top down
	the onput message o is defined from right to left, or from down to top

	each i and o is an array, aranged in a clock-wise way
"""
struct Message{_I, _O} <: AbstractMessage
	i::_I
	o::_O
end

Base.copy(x::Message) = Message(copy(x.i), copy(x.o))
Message(;i::T, o::T) where {T} = Message(i, o)


abstract type MessageNormalizationAlgorithm end
struct FixedNorm <: MessageNormalizationAlgorithm end
struct FixedSum <: MessageNormalizationAlgorithm end