package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"simulator"
	"strconv"
	"vector"
)

import "runtime/pprof"

var _ = fmt.Printf
var _ = os.Exit
var _ = pprof.StartCPUProfile

const (
	Width                  = 920
	Height                 = 720
	FPS                    = 30
	dt                     = 1 / 800.0 // s
	gravity                = 9.81 * m  // m/s^2
	m                      = 10.0
	cm                     = m / 100.0
	g                      = 1.0 // grams
	kg                     = 1000 * g
	newton                 = 1000 * g * m
	mu                     = 1.0e-3 * kg / m // Dynamic visosity, Ns/m^2
	pressureConstant       = 0.001 * m * m
	surfaceTensionConstant = 0.1

	kernelRadius          = 1.5 * cm
	particleMass          = 1.0 * g
	particleRadius        = 0.3 * cm
	collisionDampingRatio = 0.9
)

var (
	frame         = 0
	collisionMesh = createCollisionMesh()
	cpus          = runtime.NumCPU()
	rho           = 0.0 // Recalculated as needed
)

var particles ParticleList

func main() {
//	// Enable profiling
//	f, err := os.Create("profile.data")
//	if err != nil {
//		log.Fatal(err)
//	}
//	pprof.StartCPUProfile(f)
//	defer pprof.StopCPUProfile()
	particles     = NewSliceParticleList()

	// Use multiple cores.
	fmt.Printf("Using %d CPUs.\n", cpus)
	runtime.GOMAXPROCS(cpus)

	// Initialize simulator
	if err := simulator.Init("Smoothed Particle Hydrodynamics", Width, Height, FPS); err != nil {
		log.Fatalf("%v\n", err)
		return
	}
	defer simulator.Terminate()

	// Use textures for walls and ground
	simulator.TextureBoundaries("data/table.tga", "data/glass.tga")

	// Set up all the particles
	initSimulation()

	// Run the simulation main loop
	for simulator.Running() {
		simulator.WaitForNextFrame()
		simulator.Update()

		if !simulator.Paused() {
			// Update our own portion of the simulation
			updateSimulation()

			surfaceMesh := constructSurface(particles, cpus)
			simulator.AddMesh(surfaceMesh)
		}

		simulator.Draw()
	}
}

// Initialize all particles and other simulation objects.
func initSimulation() {
	// Create a grid of particles.
	//particle := simulator.LoadObjFile("data/cube.obj")["Cube"]
	particle := simulator.LoadObjFile("data/sphere.obj")["Sphere"]
	particle.Scale(particleRadius)
	spacing := 1.0
	rho = 1.0 / math.Pow(spacing, 3) * g / (cm * cm * cm) // g/cm^3
	counter := 0
	for x := -4.0; x <= 4; x += spacing {
		for y := -4.0; y <= 4; y += spacing {
			for z := 8.0; z <= 15; z += spacing {
				newParticle := new(Particle)
				newParticle.position = vector.Vector{x * cm, y * cm, z * cm}
				newParticle.velocity = vector.Zero
				newParticle.mesh = particle.Copy()
				newParticle.mesh.Translate(float32(x*cm), float32(y*cm), float32(z*cm))
				newParticle.mesh.Name += strconv.Itoa(counter)
				newParticle.density = rho
				simulator.AddMesh(newParticle.mesh)
				particles.Add(newParticle)
				counter++
			}
		}
	}

	fmt.Printf("Created %d particles.\n", counter)
}

// Create the mesh with which particles register collisions.
func createCollisionMesh() *simulator.Mesh {
	//mesh := simulator.LoadObjFile("data/openCube.obj")["OpenCube"]
	mesh := simulator.LoadObjFile("data/bowl.obj")["Bowl"]
	mesh.Scale(0.8)
	mesh.Translate(0, 0, float32(12*cm))
	mesh.ApplyTransformsToVertices()
	mesh.Name = "Boundary"
	simulator.AddMesh(mesh)
	return mesh
}

// Update our simulation
func updateSimulation() {
	frame++
	particles.ForEach(computeUpdate, cpus)
	particles.ForEach(updateParticleVelocities, cpus)
	particles.ForEach(computeCollisions, cpus)
	particles.ForEach(updateParticlePositions, cpus)
}

// Compute one step in the simulation for the particle. However, do not update
// so that other particles continue using the previous values of density, etc.
func computeUpdate(particle *Particle, cpu int) {
	// Get the neighbors and distances to them
	neighbors, distances := particles.FindNeighbors(particle, kernelRadius, cpu)

	// Initialize all accumulated quantities to zero
	density := 0.0
	colorLaplace := 0.0
	colorGradient := vector.Zero
	force := vector.Zero

	// Compute pressure due to the particle itself. This is used for computing the 
	// actual pressure by averaging this with the pressures due to other particles.
	selfPressure := pressureConstant * (particle.density - rho)

	// For all particles within the kernel radius
	for i, neighbor := range neighbors {
		// Density. Note that unlike other properties, the particle contributes
		// to its own density.  This avoids having particles with a zero density.
		density += particleMass * smoothingKernel(distances[i])

		if neighbor != particle {
			// Pressure
			neighborPressure := pressureConstant * (neighbor.density - rho)
			pressureScale := particleMass * (selfPressure + neighborPressure) / (2 * neighbor.density) * derivPressureKernel(distances[i])

			neighborDirection := neighbor.position.Subtract(particle.position)
			force.Translate(neighborDirection.Scale(pressureScale / distances[i]))

			// Surface tension (color!)
			colorScaling := particleMass / neighbor.density * derivSmoothingKernel(distances[i])
			colorLaplace += particleMass / neighbor.density * laplaceSmoothingKernel(distances[i])
			colorGradient.Translate(neighborDirection.Scale(colorScaling))

			// Viscosity
			viscosityScale := mu * particleMass / neighbor.density * laplaceViscosityKernel(distances[i])
			viscosityDirection := neighbor.velocity.Subtract(particle.velocity)
			force.Translate(viscosityDirection.Scale(viscosityScale))
		}
	}

	// Only count color when it exists Otherwise, this causes the normalization
	// to fail and then particles become NaN.
	if colorGradient.Length() > 0 {
		force.Translate(colorGradient.Normalized().Scale(colorLaplace * surfaceTensionConstant))
	}

	// Gravity, pointing down
	force.Z += -gravity

	// Set acceleration based on force, via F = ma
	particle.acceleration = force.Scale(1.0 / particleMass)

	// Set the to-be density
	particle.nextDensity = density
}

// Update particle state and display properties.
func updateParticleVelocities(particle *Particle, cpu int) {
	// Update velocity
	particle.velocity.Translate(particle.acceleration.Scale(dt))
}

// Update particle state and display properties.
func updateParticlePositions(particle *Particle, cpu int) {
	// Update velocity, position
	particle.position.Translate(particle.velocity.Scale(dt))

	// Update density
	particle.density = particle.nextDensity

	// Move mesh to particle location
	particle.mesh.MoveTo(particle.position)

	// If this particle is too damn far away, get rid of the bastard
	if particle.position.Z < 0 {
		particles.Remove(particle)
	}
}

// Update particle velocity due to collisions with the collision mesh.
func computeCollisions(particle *Particle, cpu int) {
	// Check if this particle has collided, and get the normal vector to the
	// surface it collided with if it has indeed collided with some face.
	collision, normal := CheckForCollision(particle)

	if collision {
		// Upon collision, reflect the particle from the surface using the
		// normal vector. The new velocity is given by:
		//  v <- normalized(v - 2p) * ||v|| * damping, where p = proj_n v
		proj := particle.velocity.ProjectOnto(normal)
		newVelocity := particle.velocity.Subtract(proj).Subtract(proj)
		particle.velocity = newVelocity.Normalized().Scale(particle.velocity.Length() * collisionDampingRatio)
	}
}

// Check whether a particle has undergone a collision. Return whether a
// collision has happened and a vector equal to the normal vector of the face
// with which the collision happened (if the collision did happen).
func CheckForCollision(particle *Particle) (bool, vector.Vector) {
	fullNormal := vector.Zero
	collision := false
    location := particle.position
    nextLocation := particle.position.Add(particle.velocity.Scale(dt))

	// Check each face for collision separately.
	for faceInd, face := range collisionMesh.Faces {
		// Get the normal to return in case there's a collision.
		normal := collisionMesh.FaceNormals[faceInd]

		// Check a triangle.
		if len(face) == 3 {
			p0 := collisionMesh.Vertices[face[0]]
			p1 := collisionMesh.Vertices[face[1]]
			p2 := collisionMesh.Vertices[face[2]]

			// Check that this is inside the given triangle. Move the location
			// to be relative to the vertex of the triangle, so that all the
			// vectors have the same origin.
			v1 := p1.Subtract(p0)
			v2 := p2.Subtract(p0)
			if CheckTriangleCollision(location.Subtract(p0), nextLocation.Subtract(p0), v1, v2, normal) {
                collision = true
				fullNormal = fullNormal.Add(normal)
			}
		} else if len(face) == 4 {
			// Check a quadrilateral.
			p0 := collisionMesh.Vertices[face[0]]
			p1 := collisionMesh.Vertices[face[1]]
			p2 := collisionMesh.Vertices[face[2]]
			p3 := collisionMesh.Vertices[face[3]]

			// Check the first triangle in the quadrilateral.
			v1 := p1.Subtract(p0)
			v2 := p3.Subtract(p0)
			if CheckTriangleCollision(location.Subtract(p0), nextLocation.Subtract(p0), v1, v2, normal) {
				collision = true
				fullNormal = fullNormal.Add(normal)
			}

			// Check the second triangle in the quadrilateral.
			v1 = p1.Subtract(p2)
			v2 = p3.Subtract(p2)
			if CheckTriangleCollision(location.Subtract(p2), nextLocation.Subtract(p2), v1, v2, normal) {
				collision = true
				fullNormal = fullNormal.Add(normal)
			}
		}
	}

	// Normalize the total vector if a collision happened.
	if collision {
		return true, fullNormal.Normalized()
	}

	// No collision happened with any face.
	return false, vector.Zero
}

// Check if the triangle defined by two vectors v1 and v2 and a normal vector
// norm has a collision (within a particle radius) of the given location.
func CheckTriangleCollision(location, nextLocation, v1, v2, norm vector.Vector) bool {
    // Compute the ray from the current location to the next location
    rayStart := location
    rayDirection := nextLocation.Subtract(location)
    intersectionDist := -norm.Dot(rayStart) / norm.Dot(rayDirection)

    // Check that the location is close enough to the face to intersect in the next timestep
    if intersectionDist < 0 || intersectionDist > 1 {
        return false
    }

    // Find the intersection
    intersection := rayStart.Add(rayDirection.Scale(intersectionDist))

    // Check if it is inside the triangle
    if norm.Dot(v1.CrossProduct(intersection)) < 0 {
        return false
    }
    if norm.Dot(v2.Subtract(v1).CrossProduct(intersection.Subtract(v1))) < 0 {
        return false
    }
    if norm.Dot(v2.Scale(-1).CrossProduct(intersection.Subtract(v2))) < 0 {
        return false
    }

    // Collision!
    return true
}
