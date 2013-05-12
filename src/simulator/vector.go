package simulator

import "math"

type Vector struct {
	X, Y, Z float64
}

// Compute length of a vector
func (v Vector) Length() float64 {
	return math.Hypot(v.X, math.Hypot(v.Y, v.Z))
}

// Normalize a vector
func Normalize(v Vector) Vector {
	length := v.Length()
	return Vector{v.X / length, v.Y / length, v.Z / length}
}

// Convert degrees to radians
func Radian(degree float64) float64 {
	return degree * math.Pi / 180
}

// Compute a dot product of two vectors.
func (vec Vector) Dot(other Vector) float64 {
	return vec.X * other.X + vec.Y * other.Y + vec.Z * other.Z
}

func (vec Vector) Subtract(other Vector) Vector {
    return Vector{vec.X - other.X, vec.Y - other.Y, vec.Z - other.Z}
}
