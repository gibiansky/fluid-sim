package vector

import "math"

// Vectors, with x, y, z components.
type Vector struct {
	X, Y, Z float64
}

// Constant zero vector.
var Zero = Vector{0, 0, 0}

// Compute distance to another vector.
func (this *Vector) DistanceTo(other Vector) float64 {
	return math.Hypot(this.X-other.X, math.Hypot(this.Y-other.Y, this.Z-other.Z))
}

// Compute vector length.
func (this *Vector) Length() float64 {
	return this.DistanceTo(Vector{0, 0, 0})
}

// Return a normalized version of this vector.
func (this Vector) Normalized() Vector {
	return this.Scale(1.0 / this.Length())
}

// Return a scaled version of this vector.
func (this Vector) Scale(scale float64) Vector {
	return Vector{this.X * scale, this.Y * scale, this.Z * scale}
}

// Return the vector difference.
func (vec Vector) Subtract(other Vector) Vector {
	return Vector{vec.X - other.X, vec.Y - other.Y, vec.Z - other.Z}
}

// Return the vector sum.
func (vec Vector) Add(other Vector) Vector {
	return Vector{vec.X + other.X, vec.Y + other.Y, vec.Z + other.Z}
}

// Compute the dot product of two vectors.
func (vec Vector) Dot(other Vector) float64 {
	return vec.X*other.X + vec.Y*other.Y + vec.Z*other.Z
}

// Compute a cross product of two vectors.
func (first Vector) CrossProduct(second Vector) Vector {
    normal := Vector{
        first.Y*second.Z - first.Z*second.Y,
        first.Z*second.X - first.X*second.Z,
        first.X*second.Y - first.Y*second.X,
    }
	return normal
}

// Compute the projection of this vector onto the argument vector.
func (vec Vector) ProjectOnto(other Vector) Vector {
	scale := vec.Dot(other) / math.Pow(other.Length(), 2)
	return other.Scale(scale)
}

// Shift this vector by some translation. This modifies the original vector.
func (vec *Vector) Translate(other Vector) {
	vec.X += other.X
	vec.Y += other.Y
	vec.Z += other.Z
}
