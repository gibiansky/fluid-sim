package main

import (
	"fmt"
	"math"
	"os"
	"simulator"
	"sync"
	"vector"
)

var _ = fmt.Printf
var _ = os.Exit

const (
	cubeWidth    = 0.4 * cm
	isoThreshold = 500.0
)

var (
	pointArray           [][][]float64
	minVector, maxVector vector.Vector
	mutex                sync.Mutex

	kernelRadiusVec = vector.Vector{kernelRadius, kernelRadius, kernelRadius}
)

// Construct the surface mesh for a set of particles
func constructSurface(particles ParticleList, cpus int) *simulator.Mesh {
	// Get the true bounding box
	minVector, maxVector = findBoundingBox(particles, cpus)

	// Expand bounding box to by two radii on each side to allow for padding
	minVector = minVector.Subtract(kernelRadiusVec.Scale(2))
	maxVector = maxVector.Add(kernelRadiusVec.Scale(2))

	// Allocate space for the field values
	pointArray = makePointArray(minVector, maxVector)

	// Compute field values at each vertex
	particles.ForEach(computeFieldValues, cpus)

	// Create the actual mesh
	vertices, faces := makeSurfaceMesh()

    // Optimize the mesh
    vertices, faces = optimizeMesh(vertices, faces)

	// Make the mesh object and return it
	return simulator.CreateMesh("Surface", vertices, faces)
}

// Optimize a mesh by removing vertices and joining faces
func optimizeMesh(vertices []vector.Vector, faces [][]int64) ([]vector.Vector, [][]int64) {
    fmt.Printf("Vertices: %v\nFaces: %v\n\n", len(vertices), len(faces))
    return vertices, faces
}

// Create the mesh using marching cubes. The values at each point must already
// be initialized properly.
func makeSurfaceMesh() ([]vector.Vector, [][]int64) {
	// Allocate initial arrays for vertices and faces
	vertices := make([]vector.Vector, 0, 10)
	faces := make([][]int64, 0, 10)

	// Convenience function to add a face made of three vertices
	addFace := func(v1, v2, v3 vector.Vector) {
        newVerts := 3
		numVerts := int64(len(vertices))
        n1, n2, n3 := int64(-1), int64(-1), int64(-1)

        distThresh := 0.01 * cm
        for i := len(vertices) - 1; i >= 0 && i >= int(0.8 * float32(len(vertices))); i-- {
            if v1.DistanceTo(vertices[i]) < distThresh && n1 < 0 {
                n1 = int64(i)
                newVerts --
            } else if v2.DistanceTo(vertices[i]) < distThresh && n2 < 0 {
                n2 = int64(i)
                newVerts --
            } else if v3.DistanceTo(vertices[i]) < distThresh && n3 < 0 {
                n3 = int64(i)
                newVerts --
            }
        }

        used := int64(0)
        if n1 < 0 {
            n1 = numVerts + used
            used ++
        }
        if n2 < 0 {
            n2 = numVerts + used
            used ++
        }
        if n3 < 0 {
            n3 = numVerts + used
        }

		// Create the face
		face := []int64{n1, n2, n3}

		// If there's not enough room for one more face, allocate more space
		l := len(faces)
		if l+1 > cap(faces) {
			// Allocate double what's needed, for future growth.
			newSlice := make([][]int64, (l+1)*2)
			copy(newSlice, faces)
			faces = newSlice
		}

		// Store new face
		faces = faces[0 : l+1]
		faces[l] = face

		// If there's not enough room for three more vertices, allocate more space
		l = len(vertices)
		if l+newVerts > cap(vertices) {
			// Allocate double what's needed, for future growth.
			newSlice := make([]vector.Vector, (l+newVerts)*2)
			copy(newSlice, vertices)
			vertices = newSlice
		}

		// Store new vertices
		vertices = vertices[0 : l+newVerts]
		vertices[n1] = v1
		vertices[n2] = v2
		vertices[n3] = v3
	}

	// Do the marching cubes craziness
	xpts := len(pointArray)
	ypts := len(pointArray[0])
	zpts := len(pointArray[0][0])
	for x := 0; x < xpts-1; x++ {
		for y := 0; y < ypts-1; y++ {
			for z := 0; z < zpts-1; z++ {
				// Stick grid values into an array indexed by the vertex number
				gridVal := [8]float64{
					pointArray[x][y][z],
					pointArray[x][y+1][z],
					pointArray[x+1][y+1][z],
					pointArray[x+1][y][z],
					pointArray[x][y][z+1],
					pointArray[x][y+1][z+1],
					pointArray[x+1][y+1][z+1],
					pointArray[x+1][y][z+1],
				}

				// Stick the locations of the grid points into an array (indexed similarly)
				makeCorner := func(ix, iy, iz int) vector.Vector {
					return minVector.Add(vector.Vector{float64(ix), float64(iy), float64(iz)}.Scale(cubeWidth))
				}
				gridP := [8]vector.Vector{
					makeCorner(x, y, z),
					makeCorner(x, y+1, z),
					makeCorner(x+1, y+1, z),
					makeCorner(x+1, y, z),
					makeCorner(x, y, z+1),
					makeCorner(x, y+1, z+1),
					makeCorner(x+1, y+1, z+1),
					makeCorner(x+1, y, z+1),
				}

				// Compute cube index (0 through 255)
				cubeindex := cubeIndex(gridVal)

				// If the cube is not entirely within the surface (or entirely out of the surface)
				if edgeTable[cubeindex] != 0 {
					// Compute all vertices that might be in these triangles (interpolation)
					vertlist := getVertices(cubeindex, gridP, gridVal)

					// Add the triangles
					for i := 0; triTable[cubeindex][i] != -1; i += 3 {
						v1 := vertlist[triTable[cubeindex][i]]
						v2 := vertlist[triTable[cubeindex][i+1]]
						v3 := vertlist[triTable[cubeindex][i+2]]
						addFace(v1, v2, v3)
					}
				}
			}
		}
	}

	return vertices, faces
}

// Compute the cube index
func cubeIndex(gridVal [8]float64) uint8 {
	cubeindex := 0
	if gridVal[0] < isoThreshold {
		cubeindex |= 1
	}
	if gridVal[1] < isoThreshold {
		cubeindex |= 2
	}
	if gridVal[2] < isoThreshold {
		cubeindex |= 4
	}
	if gridVal[3] < isoThreshold {
		cubeindex |= 8
	}
	if gridVal[4] < isoThreshold {
		cubeindex |= 16
	}
	if gridVal[5] < isoThreshold {
		cubeindex |= 32
	}
	if gridVal[6] < isoThreshold {
		cubeindex |= 64
	}
	if gridVal[7] < isoThreshold {
		cubeindex |= 128
	}

	return uint8(cubeindex)
}

// Get a list of vertices that are used in the triangles for this cube The
// result is a 12-length slice where the ith value is the vertex on the ith
// edge (or just the zero vector, if there is no vertex).
func getVertices(cubeindex uint8, gridP [8]vector.Vector, gridVal [8]float64) [12]vector.Vector {
	var vertlist [12]vector.Vector

	// Function to interpolate between two vertices
	vertexInterp := func(p1, p2 vector.Vector, v1, v2 float64) vector.Vector {
		if math.Abs(isoThreshold-v1) < 0.00001 {
			return p1
		}
		if math.Abs(isoThreshold-v2) < 0.00001 {
			return p2
		}
		if math.Abs(v1-v2) < 0.00001 {
			return p1
		}

		mu := (isoThreshold - v1) / (v2 - v1)

		return vector.Vector{
			p1.X + mu*(p2.X-p1.X),
			p1.Y + mu*(p2.Y-p1.Y),
			p1.Z + mu*(p2.Z-p1.Z),
		}
	}

	/* Find the vertices where the surface intersects the cube */
	if edgeTable[cubeindex]&1 != 0 {
		vertlist[0] = vertexInterp(gridP[0], gridP[1], gridVal[0], gridVal[1])
	}
	if edgeTable[cubeindex]&2 != 0 {
		vertlist[1] = vertexInterp(gridP[1], gridP[2], gridVal[1], gridVal[2])
	}
	if edgeTable[cubeindex]&4 != 0 {
		vertlist[2] = vertexInterp(gridP[2], gridP[3], gridVal[2], gridVal[3])
	}
	if edgeTable[cubeindex]&8 != 0 {
		vertlist[3] = vertexInterp(gridP[3], gridP[0], gridVal[3], gridVal[0])
	}
	if edgeTable[cubeindex]&16 != 0 {
		vertlist[4] = vertexInterp(gridP[4], gridP[5], gridVal[4], gridVal[5])
	}
	if edgeTable[cubeindex]&32 != 0 {
		vertlist[5] = vertexInterp(gridP[5], gridP[6], gridVal[5], gridVal[6])
	}
	if edgeTable[cubeindex]&64 != 0 {
		vertlist[6] = vertexInterp(gridP[6], gridP[7], gridVal[6], gridVal[7])
	}
	if edgeTable[cubeindex]&128 != 0 {
		vertlist[7] = vertexInterp(gridP[7], gridP[4], gridVal[7], gridVal[4])
	}
	if edgeTable[cubeindex]&256 != 0 {
		vertlist[8] = vertexInterp(gridP[0], gridP[4], gridVal[0], gridVal[4])
	}
	if edgeTable[cubeindex]&512 != 0 {
		vertlist[9] = vertexInterp(gridP[1], gridP[5], gridVal[1], gridVal[5])
	}
	if edgeTable[cubeindex]&1024 != 0 {
		vertlist[10] = vertexInterp(gridP[2], gridP[6], gridVal[2], gridVal[6])
	}
	if edgeTable[cubeindex]&2048 != 0 {
		vertlist[11] = vertexInterp(gridP[3], gridP[7], gridVal[3], gridVal[7])
	}

	return vertlist
}

// Find the bounding box of the particles in parallel and return the minimum and maximum corners
func findBoundingBox(particles ParticleList, cpus int) (vector.Vector, vector.Vector) {
	// Keep running minima and maxima, one per processor per coordinate
	xmin, ymin, zmin := make([]float64, cpus), make([]float64, cpus), make([]float64, cpus)
	xmax, ymax, zmax := make([]float64, cpus), make([]float64, cpus), make([]float64, cpus)

	// Initialize minima/maxima to positive/negative infinities, respectively
	for i := 0; i < cpus; i++ {
		xmin[i], ymin[i], zmin[i] = math.Inf(1), math.Inf(1), math.Inf(1)
		xmax[i], ymax[i], zmax[i] = math.Inf(-1), math.Inf(-1), math.Inf(-1)
	}

	// Compute minima and maxima in parallel for each processor
	particles.ForEach(func(particle *Particle, cpu int) {
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
	}, cpus)

	// Create the vectors representing the two opposite corners of the obunding box
	// Initialize them with values from the first processor
	minVector := vector.Vector{xmin[0], ymin[0], zmin[0]}
	maxVector := vector.Vector{xmax[0], ymax[0], zmax[0]}

	// Find minima and maxima for each coordinate across all processors
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

// Make a rectangular prism bounding box mesh (for debugging)
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

	faces[0] = []int64{0, 1, 3, 2}
	faces[1] = []int64{0, 1, 5, 4}
	faces[2] = []int64{0, 2, 6, 4}
	faces[3] = []int64{2, 3, 7, 6}
	faces[4] = []int64{1, 3, 7, 5}
	faces[5] = []int64{4, 5, 7, 6}

	return simulator.CreateMesh("BoundingBox", vertices, faces)
}

// Construct the array to hold the field values
func makePointArray(minVector, maxVector vector.Vector) [][][]float64 {
	// Compute size of array with leniency
	diffVector := maxVector.Subtract(minVector)
	cubeVector := diffVector.Scale(1 / cubeWidth).Add(vector.Vector{1, 1, 1})
	xpts, ypts, zpts := int(cubeVector.X), int(cubeVector.Y), int(cubeVector.Z)

	// Create array and all sub-arrays
	points := make([][][]float64, xpts)
	for i := 0; i < xpts; i++ {
		points[i] = make([][]float64, ypts)
		for j := 0; j < ypts; j++ {
			points[i][j] = make([]float64, zpts)
		}
	}

	return points
}

// Compute field value contributions from a particle to any points it has influence over
func computeFieldValues(particle *Particle, cpu int) {
	// Shift to have position within bounding box
	boundingPosition := particle.position.Subtract(minVector)

	// Get corners of box in which this particle can exert any influence at all (in cube coords)
	lower := boundingPosition.Subtract(kernelRadiusVec).Scale(1 / cubeWidth)
	upper := boundingPosition.Add(kernelRadiusVec).Scale(1 / cubeWidth)

	// Compute cube index positions of those corners
	xStart, yStart, zStart := math.Ceil(lower.X), math.Ceil(lower.Y), math.Ceil(lower.Z)
	xEnd, yEnd, zEnd := math.Floor(upper.X), math.Floor(upper.Y), math.Floor(upper.Z)

	// Loop over all cubes this particle might influence
	for i := xStart; i <= xEnd; i++ {
		for j := yStart; j <= yEnd; j++ {
			for k := zStart; k <= zEnd; k++ {
				// Compute distance to that cube corner
				dist := boundingPosition.DistanceTo(vector.Vector{i, j, k}.Scale(cubeWidth))

				// If the distance is within the kernel, increment by the amount of influence
				if dist < kernelRadius {
					// Do not allow multiple threads to write to the same location.
					// Otherwise, we could get two threads ignoring each others' values.
					mutex.Lock()
					pointArray[int(i)][int(j)][int(k)] += smoothingKernel(dist)
					mutex.Unlock()
				}
			}
		}
	}
}
