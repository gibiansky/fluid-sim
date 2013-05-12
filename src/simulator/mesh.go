package simulator

import (
	"fmt"
	"io/ioutil"
    "vector"
	"math"
	"strconv"
	"strings"
)

// Mesh object loaded from a *.obj file.
type Mesh struct {
	// Human-readable name for this mesh.
	Name string

	// Array of vertices used in this mesh.
	Vertices []vector.Vector

	// List of faces. Each face is a list of vertex indices.
	Faces [][]int64

	// List of face normals and vertex normals
	FaceNormals   []vector.Vector
	VertexNormals []vector.Vector

	// Transformation applied to this mesh before drawing
	dx, dy, dz, rotX, rotY, rotZ, sx, sy, sz float32
}

/*** Exposed methods ***/

// Return a copy of a mesh. The copy shares the vertex and face array.
func (original *Mesh) Copy() *Mesh {
	mesh := new(Mesh)
	mesh.sx = original.sx
	mesh.sy = original.sy
	mesh.sz = original.sz
	mesh.Name = original.Name
	mesh.Vertices = original.Vertices
	mesh.Faces = original.Faces
	mesh.FaceNormals, mesh.VertexNormals = original.FaceNormals, original.VertexNormals
	return mesh
}

// Return a deep copy of the mesh, which does not share the vertex and face
// array.  Use this only if you plan on modifying the vertices or faces;
// otherwise, use Copy().
func (original *Mesh) DeepCopy() *Mesh {
	mesh := original.Copy()

	// Copy over vertex data
	mesh.Vertices = make([]vector.Vector, len(original.Vertices))
	copy(mesh.Vertices, original.Vertices)

	// Copy over face data
	mesh.Faces = make([][]int64, len(original.Faces))
	for i := 0; i < len(original.Faces); i++ {
		mesh.Faces[i] = make([]int64, len(original.Faces[i]))
		copy(mesh.Faces[i], original.Faces[i])
	}

	return mesh
}

// Translate the mesh some amount.
func (mesh *Mesh) Translate(dx, dy, dz float32) {
	mesh.dx += dx
	mesh.dy += dy
	mesh.dz += dz
}

// Translate the mesh and apply to verts
func (mesh *Mesh) ApplyTransformsToVertices() {
    for i, _ := range mesh.Vertices {
        mesh.Vertices[i].X *= float64(mesh.sx)
        mesh.Vertices[i].Y *= float64(mesh.sy)
        mesh.Vertices[i].Z *= float64(mesh.sz)
    }

    for i, _ := range mesh.Vertices {
        mesh.Vertices[i].X += float64(mesh.dx)
        mesh.Vertices[i].Y += float64(mesh.dy)
        mesh.Vertices[i].Z += float64(mesh.dz)
    }

    mesh.sx = 1
    mesh.sy = 1
    mesh.sz = 1
    mesh.dx = 0
    mesh.dy = 0
    mesh.dz = 0
}

func (mesh *Mesh) MoveTo(loc vector.Vector) {
	mesh.dx = float32(loc.X)
	mesh.dy = float32(loc.Y)
	mesh.dz = float32(loc.Z)
}

// Rotate the mesh by some angle around each axis.  This is currently
// implemented as rotating about the x-axis, then the y-axis, then the z-axis,
// by the amounts specified.  Thus, successive Rotate calls might not actually
// rotate correctly.  
// TODO: Fix rotation of meshes.
func (mesh *Mesh) Rotate(rotX, rotY, rotZ float32) {
	mesh.rotX += rotX
	mesh.rotY += rotY
	mesh.rotZ += rotZ
}

// Translate the mesh by some factor in all directions.
func (mesh *Mesh) Scale(factor float32) {
	mesh.sx = factor
	mesh.sy = factor
	mesh.sz = factor
}

// Return the location of the mesh as a vector.
func (mesh *Mesh) Location() vector.Vector {
	return vector.Vector{float64(mesh.dx), float64(mesh.dy), float64(mesh.dz)}
}

/*** Exposed functions ***/

// Load all meshes in a *.obj file. Currently, this only supports files that
// only use simple vertices and faces, and do not have vertex normals or
// texture coordinates specified.  Returns a map from the mesh names to the
// mesh objects.
func LoadObjFile(filename string) map[string]*Mesh {
	// Read all contents of file
	contents, err := ioutil.ReadFile(filename)
	if err != nil {
		panic(fmt.Sprintf("File '%s' does not exist!", filename))
	}

	// Split file into lines.
	lines := strings.Split(string(contents), "\n")

	// Count number of objects. New objects are denoted by lines that start
	// with the character 'o'.
	numMeshes := 0
	for _, line := range lines {
		if len(line) > 0 && line[0] == 'o' {
			numMeshes++
		}
	}

	// Count vertices and faces in each mesh. Note that vertices are shared
	// between meshes; there is a single vertex array for the entire obj file.
	numVertices := 0
	numFaces := make([]int, numMeshes)
	meshIndex := 0
	for _, line := range lines {
		// Ignore the last line, which has length zero.
		if len(line) > 0 {
			// The first character determines what information this line contains.
			switch line[0] {
			case 'v':
				// Vertices are denoted by a line starting with 'v'.
				numVertices++
			case 'f':
				// Faces are denoted by a line starting with 'f'. We subtract
				// one because we've already seen one line starting with an 'o'
				// by the time we get to the first face, and arrays are
				// zero-indexed.
				numFaces[meshIndex-1]++
			case 'o':
				// New objects are denoted by a line starting with 'o'. Upon
				// encountering the new object, change which object we're
				// counting faces for.
				meshIndex++
			}
		}
	}

	// Allocate space for storing mesh data and vertex data.
	meshes := make(map[string]*Mesh)
	vertices := make([]vector.Vector, numVertices)

	var faces [][]int64
	vertexIndex := 0
	faceIndex := 0
	meshName := ""
	meshIndex = 0

	// Helper function to create a mesh object and store it.
	storeMesh := func() {
		mesh := new(Mesh)
		mesh.Name = meshName
		mesh.Vertices = vertices
		mesh.Faces = faces
		mesh.FaceNormals, mesh.VertexNormals = computeNormals(mesh)
		meshes[meshName] = mesh
		mesh.sx = 1
		mesh.sy = 1
		mesh.sz = 1
	}

	// Read vertex and face data into the meshes.
	for _, line := range lines {
		// Ignore the last line with length zero.
		if len(line) > 0 {
			// First character determines what a line contains.
			switch line[0] {
			case 'o':
				// This line contains a new object declaration, with a name.
				// If this is a new mesh, save the previous one.
				if meshIndex > 0 {
					storeMesh()
				}

				// Parse the mesh name and allocate space for the face list.
				meshName = line[2:]
				faces = make([][]int64, numFaces[meshIndex])
				faceIndex = 0
				meshIndex++

			case 'v':
				// This line contains a vertex location. Coordinates are split
				// by a space. Parse the floating point values one by one.
				coordinateStrings := strings.Split(line[2:], " ")
				vertices[vertexIndex].X, err = strconv.ParseFloat(coordinateStrings[0], 64)
				vertices[vertexIndex].Y, err = strconv.ParseFloat(coordinateStrings[1], 64)
				vertices[vertexIndex].Z, err = strconv.ParseFloat(coordinateStrings[2], 64)
				vertexIndex++
			case 'f':
				// This line contains a face specified by a space-separated list of vertex indices.
				faceStrings := strings.Split(line[2:], " ")

				// Allocate space for the list of vertex indices in this face.
				faces[faceIndex] = make([]int64, len(faceStrings))

				// Parse the index list denoting this face and store it.
				for i := 0; i < len(faceStrings); i++ {
					faces[faceIndex][i], err = strconv.ParseInt(faceStrings[i], 10, 0)

					// Switch to zero-based indexing of vertex indices. The
					// file format uses 1 to refer to the first vertex, but we
					// use 0, so subtract one.
					faces[faceIndex][i]--
				}

				faceIndex++
			}

		} else {
			// Store the last mesh in the file.
			storeMesh()
		}
	}

	return meshes
}

/*** Internal functions ***/

// Compute face and vertex normal vectors for a mesh.
func computeNormals(mesh *Mesh) ([]vector.Vector, []vector.Vector) {
	// Allocate space for face and vertex normals.
	faceNormals := make([]vector.Vector, len(mesh.Faces))
	vertexNormals := make([]vector.Vector, len(mesh.Vertices))

	// Compute normals for faces. Normals for faces are formed by the cross
	// product of their edges.
	for faceInd, face := range mesh.Faces {
		// Compute edge 1 and edge 2 as vectors between pairs of vertices
		var first vector.Vector
		first.X = mesh.Vertices[face[0]].X - mesh.Vertices[face[1]].X
		first.Y = mesh.Vertices[face[0]].Y - mesh.Vertices[face[1]].Y
		first.Z = mesh.Vertices[face[0]].Z - mesh.Vertices[face[1]].Z

		var second vector.Vector
		second.X = mesh.Vertices[face[0]].X - mesh.Vertices[face[2]].X
		second.Y = mesh.Vertices[face[0]].Y - mesh.Vertices[face[2]].Y
		second.Z = mesh.Vertices[face[0]].Z - mesh.Vertices[face[2]].Z

		// Compute Normal = CrossProduct(edge 1, edge 2)
		faceNormals[faceInd] = first.CrossProduct(second).Normalized()
	}

	// Compute normals for vertices. Normals for vertices are the sum of the
	// faces to which this vertex belonged. (Note that after taking the sum,
	// the normals must be renormalized to be unit length.)
	for vertInd, _ := range mesh.Vertices {

		// Initialize vertex normals to zero
		vertexNormals[vertInd].X, vertexNormals[vertInd].Y, vertexNormals[vertInd].Z = 0, 0, 0

		// Look at all faces to find the ones where this vertex is present.
		for faceInd, face := range mesh.Faces {
			for _, vert := range face {
				// If this vertex is present in this face, add the face normal
				// to the vertex normal. This only happens once per face max.
				if int(vert) == vertInd {
					vertexNormals[vertInd].X += faceNormals[faceInd].X
					vertexNormals[vertInd].Y += faceNormals[faceInd].Y
					vertexNormals[vertInd].Z += faceNormals[faceInd].Z
				}
			}
		}

		// Renormalize the normal to be unit length.
		normalSize := math.Hypot(vertexNormals[vertInd].X, math.Hypot(vertexNormals[vertInd].Y, vertexNormals[vertInd].Z))
		vertexNormals[vertInd].X /= normalSize
		vertexNormals[vertInd].Y /= normalSize
		vertexNormals[vertInd].Z /= normalSize
	}

	return faceNormals, vertexNormals
}
