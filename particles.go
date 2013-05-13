package main

import (
	"math"
	"simulator"
	"vector"
)

type Particle struct {
	position, velocity   vector.Vector
	mesh                 *simulator.Mesh
	acceleration         vector.Vector
	density, nextDensity float64
}

type ParticleList interface {
	Add(*Particle)
	Remove(*Particle)
	FindNeighbors(*Particle, float64, int) ([]*Particle, []float64)
	ForEach(func(*Particle, int), int)
}

type SliceParticleList struct {
	particles      []*Particle
	neighborSlices [][]*Particle
	distanceSlices [][]float64
}

func (this *SliceParticleList) Add(particle *Particle) {
	l := len(this.particles)
	if l+1 > cap(this.particles) { // reallocate
		// Allocate double what's needed, for future growth.
		newSlice := make([]*Particle, (l+1)*2)
		copy(newSlice, this.particles)
		this.particles = newSlice
	}
	this.particles = this.particles[0 : l+1]
	this.particles[l] = particle
}

// Remove a particle from the simulation
func (this *SliceParticleList) Remove(particle *Particle) {
	// Remove the particle's mesh from the collection of drawn objects
	simulator.DeleteMesh(particle.mesh.Name)

	// Remove the particle from our collection to update
	index := -1
	for i, p := range this.particles {
		if particle == p {
			index = i
			break
		}
	}
	this.particles[index] = this.particles[len(this.particles)-1]
	this.particles = this.particles[0 : len(this.particles)-1]
}

func NewSliceParticleList() *SliceParticleList {
	thing := new(SliceParticleList)
	return thing
}

func (this *SliceParticleList) FindNeighbors(particle *Particle, radius float64, proc int) ([]*Particle, []float64) {
	counter := 0
	for _, p := range this.particles {
		// Optimize a bit by checking point-wise distance before computing true distance
		delta := particle.position.Subtract(p.position)
		if math.Abs(delta.X) < radius && math.Abs(delta.Y) < radius && math.Abs(delta.Z) < radius {

			// Since we know there's a possibility of this being within the radius, compute the precise distance
			this.distanceSlices[proc][counter] = particle.position.DistanceTo(p.position)
			if this.distanceSlices[proc][counter] <= radius {
				this.neighborSlices[proc][counter] = p
				counter++
			}
		}
	}
	return this.neighborSlices[proc][0:counter], this.distanceSlices[proc][0:counter]
}

func (this *SliceParticleList) ForEach(todo func(*Particle, int), processors int) {
	// If there are not the right about of storage slices (one per processor), 
	// reallocate space to store neighbors and distances.
	if len(this.neighborSlices) != processors {
		this.neighborSlices = make([][]*Particle, processors)
		this.distanceSlices = make([][]float64, processors)

		// Allocate each slice
		for i, _ := range this.neighborSlices {
			this.neighborSlices[i] = make([]*Particle, len(this.particles))
			this.distanceSlices[i] = make([]float64, len(this.particles))
		}
	}

	// Start a new goroutine for each processor
	particlesPerProcessor := len(this.particles) / processors
	done := make(chan bool)
	for proc := 0; proc < processors; proc++ {
		start := proc * particlesPerProcessor
		end := (proc + 1) * particlesPerProcessor
		if proc == processors-1 {
			end = len(this.particles)
		}

		go func(proc int) {
			for _, particle := range this.particles[start:end] {
				todo(particle, proc)
			}

			done <- true
		}(proc)
	}

	for i := 0; i < processors; i++ {
		<-done
	}
}
