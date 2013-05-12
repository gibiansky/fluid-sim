package main

import (
	"fmt"
	"math"
	"os"
	"simulator"
	"vector"
)


var _ = fmt.Printf
var _ = os.Exit

const (
	cubeWidth = 0.5 * cm
)

var (
)


func constructSurface(particles ParticleList, cpus int) *simulator.Mesh {
	// construct bounding box
	minVector, maxVector := findBoundingBox(particles, cpus)

	return makeBoundingBoxMesh(minVector, maxVector)

	// construct mesh of cubes


	// compute field values at each vertex
}

func findBoundingBox(particles ParticleList, cpus int) (vector.Vector, vector.Vector) {
	// Keep running min/max.
	xmin, ymin, zmin := make([]float64, cpus), make([]float64, cpus), make([]float64, cpus)
	xmax, ymax, zmax := make([]float64, cpus), make([]float64, cpus), make([]float64, cpus)

	for i := 0; i < cpus; i++ {
		xmin[i], ymin[i], zmin[i] = math.Inf(1), math.Inf(1), math.Inf(1)
		xmax[i], ymax[i], zmax[i] = math.Inf(-1), math.Inf(-1), math.Inf(-1)
	}
	particles.ForEach(func (particle *Particle, cpu int) {
		if particle.position.X > xmax[cpu] {
			xmax[cpu] = particle.position.X
		}
		if particle.position.Y > ymax[cpu] {
			ymax[cpu] = particle.position.Y
		}
		if particle.position.Z > zmax[cpu] {
			zmax[cpu] = particle.position.Z
		}

		if particle.position.X < xmin[cpu] {
			xmin[cpu] = particle.position.X
		}
		if particle.position.Y < ymin[cpu] {
			ymin[cpu] = particle.position.Y
		}
		if particle.position.Z < zmin[cpu] {
			zmin[cpu] = particle.position.Z
		}
	},cpus)

	// Find the bounding box
	minVector := vector.Vector{xmin[0], ymin[0], zmin[0]}
	maxVector := vector.Vector{xmax[0], ymax[0], zmax[0]}
	for i := 1; i < cpus; i++ {
		minVector.X = math.Min(minVector.X, xmin[i])
		minVector.Y = math.Min(minVector.Y, ymin[i])
		minVector.Z = math.Min(minVector.Z, zmin[i])

		maxVector.X = math.Max(maxVector.X, xmax[i])
		maxVector.Y = math.Max(maxVector.Y, ymax[i])
		maxVector.Z = math.Max(maxVector.Z, zmax[i])
	}

	return minVector, maxVector

}

func makeBoundingBoxMesh(minVector, maxVector vector.Vector) *simulator.Mesh {
	vertices := make([]vector.Vector, 8)
	faces := make([][]int64, 6)

	vertices[0] = vector.Vector{minVector.X, minVector.Y, minVector.Z}
	vertices[1] = vector.Vector{minVector.X, minVector.Y, maxVector.Z}
	vertices[2] = vector.Vector{minVector.X, maxVector.Y, minVector.Z}
	vertices[3] = vector.Vector{minVector.X, maxVector.Y, maxVector.Z}
	vertices[4] = vector.Vector{maxVector.X, minVector.Y, minVector.Z}
	vertices[5] = vector.Vector{maxVector.X, minVector.Y, maxVector.Z}
	vertices[6] = vector.Vector{maxVector.X, maxVector.Y, minVector.Z}
	vertices[7] = vector.Vector{maxVector.X, maxVector.Y, maxVector.Z}

	faces[0] = []int64{0, 1, 3, 2};
	faces[1] = []int64{0, 1, 5, 4};
	faces[2] = []int64{0, 2, 6, 4};
	faces[3] = []int64{2, 3, 7, 6};
	faces[4] = []int64{1, 3, 7, 5};
	faces[5] = []int64{4, 5, 7, 6};

	return simulator.CreateMesh("BoundingBox", vertices, faces)
}
