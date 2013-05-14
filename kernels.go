package main

import "math"

const (
    kernelRadius2 = kernelRadius * kernelRadius
    kernelRadius3 = kernelRadius2 * kernelRadius
    kernelRadius6 = kernelRadius3 * kernelRadius3
    kernelRadius9 = kernelRadius6 * kernelRadius3

)

/** Smoothing Kernels:
  - General Kernel: used for most SPH computations
  - Pressure Kernel: used for computing pressure
  - Viscosity Kernel: used for computing acceleration due to viscosity
  **/

// General kernel
func smoothingKernel(radius float64) float64 {
	return (315 / (64 * math.Pi * kernelRadius9)) * math.Pow(kernelRadius2-(radius * radius), 3)
}

// Derivative of general kernel
func derivSmoothingKernel(radius float64) float64 {
	return -(315 * 3 / (64 * math.Pi * kernelRadius9)) * math.Pow(kernelRadius2-(radius * radius), 2) * 2 * radius
}

// Laplacian of general smoothing kernel
func laplaceSmoothingKernel(radius float64) float64 {
	return -945 / 32 / kernelRadius9 / math.Pi * (3*kernelRadius2 - 7*(radius * radius)) * (kernelRadius2 - (radius * radius))
}

// Pressure kernel
func pressureKernel(radius float64) float64 {
	return 15 / (math.Pi * kernelRadius6) * math.Pow(kernelRadius-radius, 3)
}

// Gradient of pressure kernel
func derivPressureKernel(radius float64) float64 {
	return -45 / (math.Pi * kernelRadius6) * math.Pow(kernelRadius-radius, 2)
}

// Viscosity kernel
func viscosityKernel(radius float64) float64 {
    r2 := radius * radius
	ker := -((radius * r2) / (2 * kernelRadius3)) + (r2 / (kernelRadius2)) + (kernelRadius / (2 * radius)) - 1
	return ker * 15 / 2 / kernelRadius3 / math.Pi
}

// Laplacian of viscosity kernel
func laplaceViscosityKernel(radius float64) float64 {
	ker := 6 * (kernelRadius - radius) / kernelRadius3
	return ker * 15 / 2 / kernelRadius3 / math.Pi
}
